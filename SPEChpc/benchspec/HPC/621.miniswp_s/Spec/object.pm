$benchnum  = '621';
$benchname = 'miniswp_s';
$exename   = 'sweep';
$benchlang = 'C';
@base_exe  = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;

@sources = qw( 
1_base/arguments.c
1_base/env.c
1_base/env_assert.c
1_base/env_cuda.c
1_base/env_mpi.c
1_base/pointer.c
2_sweeper_base/array_operations.c
2_sweeper_base/dimensions.c
3_sweeper/faces_kba.c
3_sweeper/quantities.c
3_sweeper/stepscheduler_kba.c
3_sweeper/sweeper.c
3_sweeper/sweeper_kernels.c
4_driver/runner.c
4_driver/sweep.c
specmpitime/specmpitime.c
          );

$need_math = 'no';

$bench_cflags = '-DSPEC -DUSE_MPI -DUSE_KBA -DUSE_ACCELDIR -I1_base/ -I2_sweeper_base/ -I3_sweeper/ -Ispecmpitime';


sub invoke {
    my ($me) = @_;
    my @rc;
    my @temp = main::read_file('control');
    my $exe = $me->exe_file;
    for my $line (@temp) {
        chomp $line;
        next if $line =~ m/^\s*(?:#|$)/;
        $line =~ /(\d) (.*)/;
        my $run = $1;
        my $arg = $2;
        $arg .= " --nthread_e $me->{'threads'}";
        push (@rc, { 'command' => $exe,
                     'args'    => [ $arg ],
                     'output'  => "$exename\_$run.out",
                     'error'   => "$exename\_$run.err",
              });
    }
    return @rc;
}

%deps = (
);


sub pre_build {
        my ($me, $dir, $setup, $target) = @_;
        my $bname = "$benchnum\.$benchname";
        my $pmodel = $me->pmodel;
        $me->pre_build_set_pmodel($dir, $setup, $target, $bname);
        return 0;
}
1;
