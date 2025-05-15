$benchnum  = '613';
$benchname = 'soma_s';
$exename   = 'soma';
$benchlang = 'C';
@base_exe  = ($exename);

$bench_cflags = "-Ispecmpitime";
$reltol = 0.006;
$abstol = 0.0000004;

@sources = qw( 
allocator.c
ana.c
autotuner.c
bond.c
cmdline.c
device.c
generate_positions.c
independent_sets.c
init.c
io.c
mc.c
mesh.c
monomer.c
mpiroutines.c
phase.c
polymer.c
rng.c
soma.c
soma_config.c
soma_util.c
test.c
specmpitime/specmpitime.c
               );

$need_math = 'yes';

$bench_flags = '';

sub invoke {
    my ($me) = @_;
    my @rc;
    my @temp = main::read_file('control');
    my $exe = $me->exe_file;
    for my $line (@temp) {
	chomp $line;
	next if $line =~ m/^\s*(?:#|$)/;
	push (@rc, { 'command' => $exe, 
		     'args'    => [ $line ], 
		     'output'  => "$exename.out",
		     'error'   => "$exename.err",
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
