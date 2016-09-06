/* Source file to swap out and test different forms of log and norm_rand.
   changing the definition of LOG can point it to a different implementation.
 */
#ifndef swappableDefinitions_h
#define swappableDefinitions_h 1

#ifndef __APPLE__
  // This function is faster on windows
  #include "libm_log.h"
  #define LOG libmlog
#else
  // The libm implementation of log on Mac OSX is plenty fast.
  #define LOG log 
#endif



#endif