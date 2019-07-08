MODULE benc
  IMPLICIT NONE

  ! ==================================================================
  ! ==                 BENCHMARK OPTIONS                            ==
  ! ==                                                              ==
  ! == IBENCH(1)  : 0  default, no special actions                  ==
  ! ==              1  skip writing the restart files               ==
  ! == IBENCH(2)  : 0  default, only print timings up to 1/1000     ==
  ! ==              1  print every timer triggered                  ==
  ! ==                                                              ==
  ! ==================================================================
  INTEGER, PARAMETER ::   nbentr=10

  INTEGER :: ibench(nbentr)
  ! ==================================================================

END MODULE benc
