$benchnum  = '635';
$benchname = 'weather_s';
$exename   = 'weather';
$benchlang = 'F';
@base_exe = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;
$bench_cflags = '-Ispecmpitime';

@sources = qw(
           miniWeather.F90
	   specmpitime/specmpitime_mod.F90
	   specmpitime/specmpitime.c
);	     

$need_math = 'no';

$bench_flags = '';
$bench_fflags = '';

sub invoke {
	my ($me) = @_;
	my @rc;
        my @temp = main::read_file('control');
	my $exe = $me->exe_file;
        my $num = 0;
        for my $line (@temp) {
            ++$num;
	    push (@rc, { 'command' => $exe, 
			'args'    => [$line], 
			'output'  => "$exename.$num.stdout.out",
			'error'   => "$exename.$num.stderr.err",
			});
        }
	return @rc;
}


%deps = (
            'miniWeather.F90' => [
                'specmpitime/specmpitime_mod.F90'
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
