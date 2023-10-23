/*
libfive: a CAD kernel for modeling with implicit functions
Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <thread>

#include "libfive/tree/tree.hpp"

namespace libfive {
class Deck; /*  Forward declaration */

class BaseEvaluator
{
public:
    // Sets the vars on the oracles in the deck.
    BaseEvaluator(std::shared_ptr<Deck> deck,
                  const std::map<Tree::Id, float>& vars);
    /*Virtual destructor allows dynamic casts to other evaluator types*/    

    virtual ~BaseEvaluator()
    {
      auto          destructionTime      = std::chrono::high_resolution_clock::now();
      auto          microseconds = std::chrono::duration_cast<std::chrono::microseconds>(destructionTime - creationTime);
      std::ofstream logfile(logTimingFileName, std::ios::app | std::ios::out);
      auto          lifetime_us = microseconds.count();

      logfile << "{" << std::endl;
      logfile << "\t\"evaluator_name\":" << name << "," << std::endl;
      logfile << "\t\"evaluator_id\":" << pointerToThis() << "," << std::endl;
      logfile << "\t\"thread_id\":" << std::this_thread::get_id() << "," << std::endl;
      logfile << "\t\"z_height\":" << std::to_string(zHeight) << "," << std::endl;

      logfile << "\t\"time\": {" << std::endl;
      logfile << "\t\t\"total_us\":" << std::to_string(lifetime_us) << ", \"total_percent\":" << std::to_string(100 * lifetime_us / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"interval_us\":" << std::to_string(timeIntervalQueries) << ", \"interval_percent\":" << std::to_string(100 * timeIntervalQueries / lifetime_us) << "," << std::endl;
      // logfile << "\t\t\"array_us\":" << std::to_string(timeArrayQueries) << ", \"array_percent\":" << std::to_string(100 * timeArrayQueries / lifetime_us) << "," << std::endl;
      //logfile << "\t\t\"value_us\":" << std::to_string(timeValueQueries) << ", \"value_percent\":" << std::to_string(100 * timeValueQueries / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"values_us\":" << std::to_string(timeValuesQueries) << ", \"values_percent\":" << std::to_string(100 * timeValuesQueries / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"features_us\":" << std::to_string(timeFeaturesQueries) << ", \"features_percent\":" << std::to_string(100 * timeFeaturesQueries / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"derivs_us\":" << std::to_string(timeDerivsQueries) << ", \"derivs_percent\":" << std::to_string(100 * timeDerivsQueries / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"set_us\":" << std::to_string(timeSetQueries) << ", \"set_percent\":" << std::to_string(100 * timeSetQueries / lifetime_us) << "," << std::endl;
      logfile << "\t\t\"isinside_us\":" << std::to_string(timeIsInsideQueries) << ", \"isinside_percent\":" << std::to_string(100 * timeIsInsideQueries / lifetime_us) << "," << std::endl;

      auto remainder = lifetime_us - timeIntervalQueries /* - timeArrayQueries- timeValueQueries -*/ - timeValuesQueries -
                       timeFeaturesQueries - timeDerivsQueries - timeSetQueries - timeIsInsideQueries;
      logfile << "\t\t\"remainder_us\":" << std::to_string(remainder) << ", \"remainder_percent\":" << std::to_string(100 * remainder / lifetime_us) << std::endl;
      
      logfile << "\t}," << std::endl;
      logfile << "\t\"call_count\": {" << std::endl;
      logfile << "\t\t\"interval\":" << std::to_string(numIntervalQueries) << "," << std::endl;
      // logfile << "\t\t\"array\":" << std::to_string(numArrayQueries) << "," << std::endl;
      logfile << "\t\t\"value\":" << std::to_string(numValueQueries) << "," << std::endl;
      logfile << "\t\t\"values\":" << std::to_string(numValuesQueries) << "," << std::endl;
      logfile << "\t\t\"features\":" << std::to_string(numFeaturesQueries) << "," << std::endl;
      logfile << "\t\t\"derivs\":" << std::to_string(numDerivsQueries) << "," << std::endl;
      logfile << "\t\t\"set\":" << std::to_string(numSetQueries) << "," << std::endl;
      logfile << "\t\t\"isinside\":" << std::to_string(numIsInsideQueries) << std::endl;
      logfile << "\t}" << std::endl;


      logfile << "}," << std::endl;
    }
    
    void setName(std::string name)
    {
      this->name = name;
    }

    void* pointerToThis()
    {
      return this;
    } 
protected:
    std::shared_ptr<Deck> deck;

    void logIntervalQuery(long long duration) {
      numIntervalQueries++;
      timeIntervalQueries += duration;
    }
    /*void logArrayQuery(long long duration) {
      numArrayQueries++;
      timeArrayQueries += duration;
    }*/
    void logValueQuery() {//long long duration) {
      numValueQueries++;
      // timeValueQueries += duration;
    }
    void logValuesQuery(long long duration) {
      numValuesQueries++;
      timeValuesQueries += duration;
    }
    void logFeaturesQuery() {
      numFeaturesQueries++;
    }
    
    void logFeaturesQueryDuration(long long duration) {
      timeFeaturesQueries += duration;
    }

    void logDerivsQuery(long long duration) {
      numDerivsQueries++;
      timeDerivsQueries += duration;
    }
    void logSetQuery(long long duration) {
      numSetQueries++;
      timeSetQueries += duration;
    }
    void logIsInsideQuery(long long duration) {
      numIsInsideQueries++;
      timeIsInsideQueries += duration;
    }

  float zHeight;

private:
  
  std::chrono::time_point<std::chrono::high_resolution_clock> creationTime;

  size_t numIntervalQueries = 0;
  //size_t numArrayQueries    = 0;
  size_t numValueQueries    = 0;
  size_t numValuesQueries   = 0;
  size_t numFeaturesQueries = 0;
  size_t numDerivsQueries   = 0;
  size_t numSetQueries      = 0;
  size_t numIsInsideQueries = 0;

  long long timeIntervalQueries = 0;
  //long long timeArrayQueries    = 0;
  //long long timeValueQueries    = 0;
  long long timeValuesQueries   = 0;
  long long timeFeaturesQueries = 0;
  long long timeDerivsQueries   = 0;
  long long timeSetQueries      = 0;
  long long timeIsInsideQueries = 0;

  std::string name;
  std::string logTimingFileName  = "c:/ntop_log/slice_timing_libfive.json";
};

}   // namespace libfive
