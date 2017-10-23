        logical  convertdouble
        parameter (convertdouble = .false.)
        character*11 compiletime
        parameter (compiletime='22 Aug 2017')
        character*5 npbversion
        parameter (npbversion='3.3.1')
        character*6 cs1
        parameter (cs1='mpif90')
        character*9 cs2
        parameter (cs2='$(MPIF77)')
        character*23 cs3
        parameter (cs3='#-L/usr/local/lib -lmpi')
        character*21 cs4
        parameter (cs4='#-I/usr/local/include')
        character*4 cs5
        parameter (cs5='-O2 ')
        character*3 cs6
        parameter (cs6='-O2')
        character*6 cs7
        parameter (cs7='randi8')
