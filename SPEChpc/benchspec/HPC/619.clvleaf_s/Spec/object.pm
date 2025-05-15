$benchnum  = '619';
$benchname = 'clvleaf_s';
$exename   = 'clvleaf';
$benchlang = 'F';
@base_exe  = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;

sub invoke {
    my ($me) = @_;
    my $name;
    my @rc;

    my $exe = $me->exe_file;
    for ($me->input_files_base) {
        if (($name) = m/(.*).in$/) {
             push (@rc, { 'command' => $exe, 
                         'args'    => [ ], 
#                         'output'  => "$name.out",
#                        'output'  => "energy.txt",
                         'error'   => "$name.err",
                        });
        }
    }
    return @rc;
}



#sub compare_commands {
#    my ($me) = @_; 
#    my @rc;
#    my $clover_validate_exe = "grep 'step:' clover.out";
#
#    push @rc, (
#        {   
#            'command' => $clover_validate_exe,
#            'args'    => [],
#            'output'  => 'energy.txt',
#            'error'   => "energy.err",
#        }   
#    );  
#
#    return @rc;
#}



@sources = qw( 
accelerate.F90
accelerate_kernel.F90
advec_cell_driver.F90
advec_cell_kernel.F90
advec_mom_driver.F90
advec_mom_kernel.F90
advection.F90
build_field.F90
calc_dt.F90
calc_dt_kernel.F90
clover.F90
clover_leaf.F90
data.F90
definitions.F90
field_summary.F90
field_summary_kernel.F90
flux_calc.F90
flux_calc_kernel.F90
generate_chunk.F90
generate_chunk_kernel.F90
hydro.F90
ideal_gas.F90
ideal_gas_kernel.F90
initialise_chunk.F90
initialise_chunk_kernel.F90
initialise.F90
pack_kernel.F90
parse.F90
PdV.F90
PdV_kernel.F90
read_input.F90
report.F90
reset_field.F90
reset_field_kernel.F90
revert.F90
revert_kernel.F90
start.F90
timer.F90
timestep.F90
update_halo.F90
update_halo_kernel.F90
update_tile_halo.F90
update_tile_halo_kernel.F90
visit.F90
viscosity.F90
viscosity_kernel.F90
specmpitime/specmpitime_mod.F90
specmpitime/specmpitime.c
               );

$need_math = 'no';

$bench_flags = '-Ispecmpitime';


%deps = (
          'PdV.F90' => [
                         'ideal_gas.F90',	# module
                         'update_halo.F90',	# module
                         'revert.F90',	# module
                         'PdV_kernel.F90',	# module
                         'report.F90',	# module
                         'clover.F90',	# module
                       ],
          'accelerate.F90' => [
                                'accelerate_kernel.F90',	# module
                                'clover.F90',	# module
                              ],
          'advec_cell_driver.F90' => [
                                       'advec_cell_kernel.F90',	# module
                                       'clover.F90',	# module
                                     ],
          'advec_mom_driver.F90' => [
                                      'advec_mom_kernel.F90',	# module
                                      'clover.F90',	# module
                                    ],
          'advection.F90' => [
                               'update_halo.F90',	# module
                               'advec_mom_driver.F90',	# module
                               'advec_cell_driver.F90',	# module
                               'clover.F90',	# module
                             ],
          'build_field.F90' => [
                                 'clover.F90',	# module
                               ],
          'calc_dt.F90' => [
                             'calc_dt_kernel.F90',	# module
                             'clover.F90',	# module
                           ],
          'clover.F90' => [
                            'pack_kernel.F90',	# module
                            'definitions.F90',	# module
                            'data.F90',	# module
                          ],
          'clover_leaf.F90' => [
                                 'clover.F90',	# module
			         'specmpitime/specmpitime_mod.F90'
                               ],
          'definitions.F90' => [
                                 'data.F90',	# module
                               ],
          'field_summary.F90' => [
                                   'field_summary_kernel.F90',	# module
                                   'ideal_gas.F90',	# module
                                   'clover.F90',	# module
                                 ],
          'flux_calc.F90' => [
                               'flux_calc_kernel.F90',	# module
                               'clover.F90',	# module
                             ],
          'generate_chunk.F90' => [
                                    'generate_chunk_kernel.F90',	# module
                                    'clover.F90',	# module
                                  ],
          'hydro.F90' => [
                           'reset_field.F90',	# module
                           'advection.F90',	# module
                           'flux_calc.F90',	# module
                           'accelerate.F90',	# module
                           'PdV.F90',	# module
                           'viscosity.F90',	# module
                           'timestep.F90',	# module
                           'clover.F90',	# module
                         ],
          'ideal_gas.F90' => [
                               'ideal_gas_kernel.F90',	# module
                               'clover.F90',	# module
                             ],
          'initialise.F90' => [
                                'report.F90',	# module
                                'parse.F90',	# module
                                'clover.F90',	# module
                              ],
          'initialise_chunk.F90' => [
                                      'initialise_chunk_kernel.F90',	# module
                                      'clover.F90',	# module
                                    ],
                    'parse.F90' => [
                           'clover.F90',	# module
                           'report.F90',	# module
                           'data.F90',	# module
                         ],
          'read_input.F90' => [
                                'report.F90',	# module
                                'parse.F90',	# module
                                'clover.F90',	# module
                              ],
          'report.F90' => [
                            'clover.F90',	# module
                            'data.F90',	# module
                          ],
          'reset_field.F90' => [
                                 'reset_field_kernel.F90',	# module
                                 'clover.F90',	# module
                               ],
          'revert.F90' => [
                            'revert_kernel.F90',	# module
                            'clover.F90',	# module
                          ],
          'start.F90' => [
                           'ideal_gas.F90',	# module
                           'update_halo.F90',	# module
                           'parse.F90',	# module
                           'clover.F90',	# module
                         ],
          'timestep.F90' => [
                              'definitions.F90',	# module
                              'ideal_gas.F90',	# module
                              'calc_dt.F90',	# module
                              'viscosity.F90',	# module
                              'update_halo.F90',	# module
                              'report.F90',	# module
                              'clover.F90',	# module
                            ],
          'update_halo.F90' => [
                                 'update_halo_kernel.F90',	# module
                                 'update_tile_halo.F90',	# module
                                 'clover.F90',	# module
                               ],
          'update_halo_kernel.F90' => [
                                        'data.F90',	# module
                                      ],
          'update_tile_halo.F90' => [
                                      'update_tile_halo_kernel.F90',	# module
                                      'clover.F90',	# module
                                    ],
          'update_tile_halo_kernel.F90' => [
                                             'data.F90',	# module
                                           ],
          'viscosity.F90' => [
                               'viscosity_kernel.F90',	# module
                               'clover.F90',	# module
                             ],
          'visit.F90' => [
                           'ideal_gas.F90',	# module
                           'viscosity.F90',	# module
                           'update_halo.F90',	# module
                           'clover.F90',	# module
                         ],

        );

sub pre_build {
        my ($me, $dir, $setup, $target) = @_;
        my $bname = "$benchnum\.$benchname";
        my $pmodel = $me->pmodel;
        $me->pre_build_set_pmodel($dir, $setup, $target, $bname);
        return 0;
}

1;
