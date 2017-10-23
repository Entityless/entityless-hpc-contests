---------------------------------------------------------------
    Ocean Model Typical Calculation Sample Code for PAC2016
---------------------------------------------------------------

I. File list:
    1) ocean.f90         
        Main program. !! Not allowed to modify !!. 
        Won't be changed during the contest except necessary bug fixing.
    2) mod_data.f90    
        Auxilary data processing module. !! Not allowed to modify !!. 
        Parameters and the implementation might be replaced by Organizing Commitee during the Final phase.
        Data structure and interface won't be change.
    3) mpdata_adiff.f90
        Implementation of calculation.
        Could be modified to any style source code without changing the calculation correctness.
        Extra source files could be added here to implemented the calculation.
    4) makefile
        Program compiling and linking configuration file.
        Could be modified for optimization.
        !! Must use 'ocean.f90' as the main program !!.

II. Performance evaluation:
    Run 'make' to build the main program 'ocean', then execute the binary. The output information should be like this:
    --------------------------------------------------
     [INIT] Initializing data...
     [INIT] Done!

     [RUN]  Start calculation...
          Calculation time:1259.364(s)
     [RUN]  Done!

     [VERIFY] Start verification.
     Would be implemented later.
     [VERIFY] Done!
    --------------------------------------------------
    The time seconds in 'Calculation time' is used for performance evaluation. Lower is better.

III. Result correctness verification:
    The subroutine 'verify_data()' in mod_data.f90 would do result check. It would be provided at the site of Final together with the parameters used in the Final.
    The results calculated by original sample code should be used as the reference.
