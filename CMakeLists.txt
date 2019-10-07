cmake_minimum_required(VERSION 3.10)
project(src C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(src
        adimvalue.c
        aniso.c
        axilib.c
        commbound.c
        convert.c
        Dsys.c
        DustUpdate.c
        fftw_dummy.c
        fftw_dummy.h
        fondam.h
        Force.c
        fpe.c
        global.h
        global_ex.h
        Init.c
        Interpret.c
        LowTasks.c
        main.c
        merge.c
        mp.h
        mpi_dummy.c
        mpi_dummy.h
        Output.c
        param.h
        param_noex.h
        Pframeforce.c
        Planet.c
        proto.h
        Psys.c
        rebin.c
        RungeKunta.c
        sg.h
        sgacc.c
        sgdens.c
        sginit.c
        sgkernel.c
        sgmain.c
        sgproto.h
        sgsysinit.c
        sgupdate.c
        sgvelinit.c
        sgzero.c
        SideEuler.c
        SourceEuler.c
        split.c
        Stockholm.c
        Theo.c
        TransportEuler.c
        turb.c
        types.h
        var.c
        Viscosity.c)
