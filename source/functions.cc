#include "functions.h"

double calc_time_interval(const double & space,
                         const double & vel,
                         const unsigned int & p_degree)
                
  {
    // CFL factor 0.5 (was 0.4) for slightly larger dt; increase to 0.55–0.6 if stable
    double time_int =  space / ((p_degree +1)* vel);
    rounding_digits(time_int, 4);
    return time_int;
  }

void rounding_digits(double & value, const int & digits, bool roundown_flag) 
{
    double abs_value = std::abs(value);
    double factor;
    if (abs_value == 0.0)
    {
        return;
    }
// 1. Find the exponent of the most significant digit.
    if (abs_value < 1.0)
        {
         double exponent = std::floor(std::log10(abs_value));    
         // 2. Determine the scaling factor to move the desired rounding place to the units digit.
         factor = std::pow(10.0, exponent - digits + 1);
        }
    else // abs_value >= 1.0
        {
         // 'precision' is used as the number of decimal places (e.g., 3).
        factor = std::pow(10.0, -digits);
        }
    
    double scaled_value = abs_value / factor;
     if (roundown_flag)
        value = std::floor(scaled_value) * factor;
    else
        value = std::ceil(scaled_value / factor) * factor;

}

void echo_parameters(const Parameters &params)
{
    std::cout << "velocities " << params.vp << ",  " << params.vs
                << std::endl;

    std::cout << "frequency " << params.fo 
                << std::endl;

}