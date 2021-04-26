set (amos_sources
        amos/cacai.f
        amos/cacon.f
        amos/cairy.f
        amos/casyi.f
        amos/cbesi.f
        amos/cbesk.f
        amos/cbinu.f
        amos/cbknu.f
        amos/cbuni.f
        amos/cbunk.f
        amos/ckscl.f
        amos/cmlri.f
        amos/crati.f
        amos/cs1s2.f
        amos/cseri.f
        amos/cshch.f
        amos/cuchk.f
        amos/cunhj.f
        amos/cuni1.f
        amos/cuni2.f
        amos/cunik.f
        amos/cunk1.f
        amos/cunk2.f
        amos/cuoik.f
        amos/cwrsk.f
        amos/gamln.f
        amos/i1mach.f
        amos/r1mach.f
        amos/xerror.f
)

set_source_files_properties(${amos_sources} PROPERTIES LANGUAGE Fortran)
add_library(${AMOSLIB} ${amos_sources})
