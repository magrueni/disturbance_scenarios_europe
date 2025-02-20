/********************************************************************************************
**    iLand - an individual based forest landscape and disturbance model
**    https://iland-model.org
**    Copyright (C) 2009-  Werner Rammer, Rupert Seidl
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**    This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
********************************************************************************************/

#include "bbgenerations.h"

#include <cmath>

#include <algorithm>

#define M_PI 3.141592653589793238462

/** @class BBGenerations
    @ingroup beetlemodule
    BBGenerations calculates potential bark beetle generations based on climate data (including bark temperature).
  */
BBGenerations::BBGenerations()
{
    mUseAirTempForGenerations = false; // this is the default
}


void BBGenerations::setLatitude(double lat_degree)
{
  // based on iLand: climate.cpp
  double latitude_rad = lat_degree * M_PI / 180.;
  int day_max_len;
  int day_10_5hrs, day_14_5hrs;
  if (latitude_rad>0)
    day_max_len = 182-10; // 21st of June
  else
    day_max_len = 365-10; //southern hemisphere
  
  // calculate length of day using  the approximation formulae of: http://herbert.gandraxa.com/length_of_day.aspx
  const double j = M_PI / 182.625;
  const double ecliptic = 23.439*M_PI/180.; // to rad
  double m;
  double daylength[366];
  for (int day=0;day<366;day++) {
    m = 1. - tan(latitude_rad)*tan(ecliptic*cos(j*(day+10))); // day=0: winter solstice => subtract 10 days
    m = std::min( std::max(m, 0.), 2.);
    daylength[day] = acos(1-m)/M_PI * 24.; // result in hours [0..24]
  }
  day_10_5hrs = 0;
  for (int day=day_max_len;day<366;day++) {
    if (daylength[day]<10.5) {
      day_10_5hrs = day;
      break;
    }
  }
  day_14_5hrs = 0;
  for (int day=day_max_len;day<366;day++) {
    log( std::to_string(day) + " - " + std::to_string(daylength[day]));
    if (daylength[day]<14.5) {
      day_14_5hrs = day;
      break;
    }
  }
  mDay_10_5hrs = day_10_5hrs;
  mDay_14_5hrs = day_14_5hrs;
  
}

/**
 * @brief BBGenerations::calculateGenerations
 * @param ru
 * @return the number of filial generation (i.e, main generations) + 0.5 if a sister brood develops also for the last generation
 */
double BBGenerations::calculateGenerations()
{
    calculateBarkTemperature();

    // start at the 1. of April, and wait for 140.3 degree days (with a threhsold of 8.3 degrees)

    const SDay *clim = &mClimate[90]; // April 1st=day 91
    const SDay *last_day = &mClimate[303]; // October 31=day of year: 304
    const SDay *day_too_short = &mClimate[mDay_14_5hrs]; // day that is shorter than 14.5 hours of daylight

    double dd=0.;

    while (dd<140.3 && clim<last_day) {
    dd+=std::max(clim->tmax -8.3, 0.);
        ++clim;
    }
    // now wait for a decent warm day with tmax > 16.5 degrees
    while (clim->tmax<=16.5 && clim<last_day)
        ++clim;

    mGenerations.clear();
    // start with the first generation....
    BBGeneration bbgen;
    bbgen.start_day = clim - &mClimate[0]; // which day of year is clim?
    bbgen.is_sister_brood = false;
    bbgen.gen = 1; // the first generation
    mGenerations.push_back(bbgen);

    // process the generations
    size_t i=0;
    while (i<mGenerations.size()) {
        BBGeneration bb = mGenerations[i]; // deliberately make a copy and not a ref
        const SDay *c = &mClimate[bb.start_day];
        int doy = bb.start_day;
        double base_temp = mEffectiveBarkTemp[doy];

        double t_sum=0.;
        bool added_sister_brood = false;
        while (c < last_day) {
            t_sum = (mEffectiveBarkTemp[doy]-base_temp) / 557.;
            if (t_sum>=1.) {
                if (c<day_too_short) {
                    // start a new parental generation (a full cycle), only if the current
                    // generation is also a filial generation.
                    if (bb.is_sister_brood)
                        mGenerations.push_back(BBGeneration(doy, true, bb.gen)); // sisterbrood: (true), keep counter of filial generation
                    else
                        mGenerations.push_back(BBGeneration(doy, false, bb.gen+1)); // new filial brood (false), start a new gen

                }
                break;
            } else if (t_sum>0.5 && !added_sister_brood) {
                // start a sister brood, *if* the maximum air temperature is high enough, and if the
                // length of the day > 14.5 hours
                if (c->tmax>16.5 && c<day_too_short) {
                    mGenerations.push_back(BBGeneration(doy, true, bb.gen)); // add a sister brood generation (true), keep gen. of originating brood
                    added_sister_brood = true;
                }
            }
            ++c; ++doy;
        }
        mGenerations[i].value = std::min(t_sum, 1.);
        ++i; // now process the next generation
    }


    mNSisterBroods = 0; mNFilialBroods = 0;

    for (i=0;i<mGenerations.size();++i) {
        if (mGenerations[i].value>0.6 && mGenerations[i].is_sister_brood==true)
            mNSisterBroods = std::max( mGenerations[i].gen, mNSisterBroods );

        if (mGenerations[i].value>0.6 && mGenerations[i].is_sister_brood==false)
            mNFilialBroods = mGenerations[i].gen;
    }

    // return the number of filial broods, and increase by 0.5 if a sister brood of the last filial generation has successfully developed.
    return mNFilialBroods + ( hasSisterBrood() ? 0.5 : 0.);

//    // now accumulate the generations
//    double filial_generations = 0.;
//    double sister_generations = 0.;
//    // accumulate possible number of offspring
//    const double offspring_factor = 25.; // assuming sex-ratio of 1:1 and 50 offspring per female (see p. 59)
//    int n=0;
//    int total_offspring = 0;
//    for (i=0;i<mGenerations.count();++i) {
//        if (mGenerations[i].value>0.6)  {
//            ++n;
//            total_offspring+= pow(offspring_factor, mGenerations[i].gen-1)*2.*offspring_factor;
//        }
//        if (mGenerations[i].value>0.6 && mGenerations[i].is_sister_brood==true)
//            sister_generations+=mGenerations[i].value;
//        if (mGenerations[i].value>0.6 && mGenerations[i].is_sister_brood==false)
//            filial_generations+=mGenerations[i].value;

//    }
//    if (logLevelDebug())
//        qDebug() << "rid" <<ru->id() << "filial/sister:" << filial_generations << sister_generations << "offspring:" << total_offspring << "started generations:" << n;
//    mNSisterBroods = sister_generations;
//    mNFilialBroods = filial_generations;

//    return filial_generations + (sister_generations>0?0.5:0);
}


/**
 * Calculate the bark temperatures for this year and a given resource unit.
 * Input: climate data (tmax (C), tmean (C), radiation (MJ/m2))
 * the LAI to estimate the radiation on the ground (Wh/m2)
 * Output: calculates for each day of the year the "effective" bark-temperature and saves a cumulative sum
 * Source: Schopf et al 2004: Risikoabschaetzung von Borkenkaefermassenkalamitaeten im Nationalpark Kalkalpen
 */
void BBGenerations::calculateBarkTemperature()
{
    // estimate the fraction of light on the ground (multiplier from 0..1)
    const double k = 0.5; // constant for the beer lambert function
    double ground_light_fraction = exp(-k * mLAI );

    mFrostDaysEarly=0;
    mFrostDaysLate=0;


    for (int i=0;i<365;++i) {
        const SDay *clim = &mClimate[i];

        // radiation: MJ/m2/day -> the regression uses Wh/m2/day -> conversion-factor: 1/0.0036
        double rad_wh = clim->radiation*ground_light_fraction/0.0036;
        // calc. maximum bark temperature
        double bt_max=1.656 + 0.002955*rad_wh + 0.534*clim->tmax + 0.01884 * clim->tmax*clim->tmax;
        double diff_bt=0.;


        if (bt_max>=30.4)
            diff_bt=std::max(-310.667 + 9.603 * bt_max, 0.); // degree * hours

        // mean bark temperature
        double bt_mean=-0.173+0.0008518*rad_wh+1.054* (clim->tmax + clim->tmin)/2.;

        // effective:
        double bt_sum = std::max(bt_mean - 8.3, 0.)*24.; // degree * hours

        // corrected for days with very high bark temperature (>30.4 degrees)
        double bt_sum_eff = (bt_sum - diff_bt) / 24.; // degree * days


        mEffectiveBarkTemp[i] = (i>0? mEffectiveBarkTemp[i-1] : 0.) + bt_sum_eff;

        // frost days:
        double min_temp = clim->tmin;


        if (min_temp < -15.) {
            // note: fixed now to northern latitudes...
            if (i < 180)
                mFrostDaysEarly++;
            else
                mFrostDaysLate++;
        }

    }

}
