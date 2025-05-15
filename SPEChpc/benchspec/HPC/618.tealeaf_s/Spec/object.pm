$benchnum  = '618';
$benchname = 'tealeaf_s';
$exename   = 'tealeaf';
$benchlang = 'C';
@base_exe  = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;

$need_math = 'no';

$bench_cflags = '-I2d/drivers/ -I2d/c_kernels -Ispecmpitime';

sub invoke {
    my ($me) = @_;
    my $name;
    my @rc;

    my $exe = $me->exe_file;
    for ($me->input_files_base) {
        if (($name) = m/(.*).in$/) {
#           push (@rc, { 'command' =>run, 
             push (@rc, { 'command' => $exe, 
                         'args'    => [ ], 
                         'output'  => "output",
#                        'output'  => "energy.txt",
                         'error'   => "$name.err",
                        });
        }
    }
    return @rc;
}



@sources = qw( 
2d/chunk.c
2d/comms.c
2d/diffuse.c
2d/initialise.c
2d/main.c
2d/parse_config.c
2d/profiler.c
2d/settings.c
2d/shared.c
2d/drivers/cg_driver.c
2d/drivers/cheby_driver.c
2d/drivers/eigenvalue_driver.c
2d/drivers/field_summary_driver.c
2d/drivers/halo_update_driver.c
2d/drivers/jacobi_driver.c
2d/drivers/kernel_initialise_driver.c
2d/drivers/ppcg_driver.c
2d/drivers/remote_halo_driver.c
2d/drivers/set_chunk_data_driver.c
2d/drivers/set_chunk_state_driver.c
2d/drivers/solve_finished_driver.c
2d/drivers/store_energy_driver.c
2d/c_kernels/cg.c
2d/c_kernels/cheby.c
2d/c_kernels/field_summary.c
2d/c_kernels/jacobi.c
2d/c_kernels/kernel_initialise.c
2d/c_kernels/kernel_interface.c
2d/c_kernels/local_halos.c
2d/c_kernels/pack_halos.c
2d/c_kernels/ppcg.c
2d/c_kernels/set_chunk_data.c
2d/c_kernels/set_chunk_state.c
2d/c_kernels/solver_methods.c
2d/c_kernels/store_energy.c
2d/c_kernels/diffuse_overload.c
specmpitime/specmpitime.c
          );



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
