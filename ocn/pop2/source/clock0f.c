/*  CVS: $Id: clock0f.c,v 1.4 2003/08/12 09:06:39 lhl Exp $
*/
/*
C#######################################################################
C PSTSWM Version 4.0 (12/1/94)                                         #
C  A message-passing benchmark code and parallel algorithm testbed     #
C  that solves the nonlinear shallow water equations using the spectral#
C  transform method.                                                   #
C Written by:                                                          #
C  Patrick Worley of Oak Ridge National Laboratory                     #
C  Ian Foster of Argonne National Laboratory                           #
C Based on the sequential code STSWM 2.0 by James Hack and Ruediger    #
C  Jakob of the National Center for Atmospheric Research.              #
C Research and development funded by the Computer Hardware, Advanced   #
C  Mathematics, and Model Physics (CHAMMP) program of the U.S.         #
C  Department of Energy.                                               # 
C                                                                      #
C Questions and comments should be directed to worley@msr.epm.ornl.gov #
C Please notify and acknowledge the authors in any research or         #
C publications utilizing PSTSWM or any part of the code.               #
C                                                                      #
C NOTICE: Neither the institutions nor the authors make any            #
C representations about the suitability of this software for any       #
C purpose. This software is provided "as is", without express or       #
C implied warranty.                                                    #
C#######################################################################
*/
#include <sys/time.h>

double clock0f()
{
  double t;
  struct timeval tp;

  gettimeofday(&tp, 0);
  t = tp.tv_sec + 0.000001*tp.tv_usec;
  return(t);

}

double CLOCK0F()
{
  double t;
  struct timeval tp;

  gettimeofday(&tp, 0);
  t = tp.tv_sec + 0.000001*tp.tv_usec;
  return(t);

}

double clock0f_()
{
  double t;
  struct timeval tp;

  gettimeofday(&tp, 0);
  t = tp.tv_sec + 0.000001*tp.tv_usec;
  return(t);

}

double _clock0f()
{
  double t;
  struct timeval tp;

  gettimeofday(&tp, 0);
  t = tp.tv_sec + 0.000001*tp.tv_usec;
  return(t);

}
