$benchnum  = '634';
$benchname = 'hpgmgfv_s';
$exename   = 'hpgmgfv';
$benchlang = 'C';
@base_exe  = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;

$bench_flags = "-DUSE_BICGSTAB=1 -DUSE_SUBCOMM=1 -DUSE_FCYCLES=1 -DUSE_GSRB=1 -DBLOCKCOPY_TILE_I=32 -DBLOCKCOPY_TILE_J=4 -DBLOCKCOPY_TILE_K=16 -DBOUNDARY_TILE_I=64 -DBOUNDARY_TILE_J=16 -DBOUNDARY_TILE_K=16 -DUSE_MPI=1 -Ispecmpitime";


@sources = qw( 
hpgmg-fv.c
level.c
mg.c
operators.fv4.c
solvers.c
timers.c
directives-wrapper.c
cuda-wrapper.c
offload-fns.c
specmpitime/specmpitime.c
               ) ;

$need_math = 'yes';

sub invoke {
    my ($me) = @_;
    my @rc;

    my @runs = grep { !m/^#/ } main::read_file('control');
    my $exe = $me->exe_file;

    foreach my $run (@runs) {
        my (@args) = split(/\s+/, $run);
        push (@rc, { 'command' => $exe,
                     'args'    => [ @args ],
                     'output'  => "$exename.out",
                     'error'   => "$exename.err",
                   } );
    }
    return @rc;
}


sub pre_build {
        my ($me, $dir, $setup, $target) = @_;
        my $bname = "$benchnum\.$benchname";
        my $pmodel = $me->pmodel;
        $me->pre_build_set_pmodel($dir, $setup, $target, $bname);
        return 0;
}
1;
