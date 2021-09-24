// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_TIMER_HPP
#define MCL_TIMER_HPP 1

#include <chrono> 

namespace mcl
{

class Timer
{
protected:
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> TimeType;
    typedef std::chrono::high_resolution_clock ClockType;
    bool running;
    TimeType start_time, stop_time;

public:
    Timer() : running(true), start_time(ClockType::now()) {}

    void start()
    {
        start_time = ClockType::now();
        running = true;
    }

    void stop()
    {
        stop_time = ClockType::now();
        running = false;
    }

    double get_ms()
    {
        using namespace std::chrono;
        long duration = 0;
        if (!running) { duration = duration_cast<microseconds>(stop_time-start_time).count(); }
        else { duration = duration_cast<microseconds>(ClockType::now()-start_time).count(); }
        return double(duration) * 0.001;
    }
};

} // ns mcl

#endif