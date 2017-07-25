!**********************************************
!*aggl2d_main.f90                             *
!*                                            *
!*Module definition file for program          *
!*                                            *
!*   aggl2d v 2.0                             *
!*                                            *
!*                                            *
!*Made by                                     *
!*Kaare A Sorensen,                           *
!*20.08.98-                                   *
!**********************************************
!*Description:                                *
!* This file includes procedures for file     *
!* communication, construction of side based  *
!* datastructure from elements, coloring for  *
!* for vectorization and construction of      *
!* coarser grids from the given original one  *
!* by agglomeration - for use in multigrid    *
!**********************************************
program Agg
 use Agglomerator2d
 call main()
 stop
end
