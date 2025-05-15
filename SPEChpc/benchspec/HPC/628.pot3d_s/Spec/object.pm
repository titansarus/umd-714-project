$benchnum  = '628';
$benchname = 'pot3d_s';
$exename   = 'pot3d';
$benchlang = 'F';
@base_exe  = ($exename);

$reltol = 0.006;
$abstol = 0.0000004;

$bench_fppflags  = '-Ihdf5  -Ihdf5/include -DNDEBUG ';
$bench_flags = '-Ihdf5/include -DNDEBUG ';
$bench_cflags = '-Ihdf5/include -Ilibz/include -DNDEBUG -Ispecmpitime ';

@sources = qw( 
		number_types.F90
		pot3d.F90
		zm_parse.F90
		zm_parse_modules.F90
		zm_sds.F90
		zm_sds_modules.F90
		specmpitime/specmpitime_mod.F90
		specmpitime/specmpitime.c
		hdf5/H5.c
		hdf5/H5A.c
		hdf5/H5AC.c
		hdf5/H5ACdbg.c
		hdf5/H5ACmpio.c
		hdf5/H5ACproxy_entry.c
		hdf5/H5Abtree2.c
		hdf5/H5Adense.c
		hdf5/H5Adeprec.c
		hdf5/H5Af.c
		hdf5/H5Aint.c
		hdf5/H5Atest.c
		hdf5/H5B.c
		hdf5/H5B2.c
		hdf5/H5B2cache.c
		hdf5/H5B2dbg.c
		hdf5/H5B2hdr.c
		hdf5/H5B2int.c
		hdf5/H5B2internal.c
		hdf5/H5B2leaf.c
		hdf5/H5B2stat.c
		hdf5/H5B2test.c
		hdf5/H5Bcache.c
		hdf5/H5Bdbg.c
		hdf5/H5C.c
		hdf5/H5CS.c
		hdf5/H5Cdbg.c
		hdf5/H5Cepoch.c
		hdf5/H5Cimage.c
		hdf5/H5Clog.c
		hdf5/H5Clog_json.c
		hdf5/H5Clog_trace.c
		hdf5/H5Cmpio.c
		hdf5/H5Cprefetched.c
		hdf5/H5Cquery.c
		hdf5/H5Ctag.c
		hdf5/H5Ctest.c
		hdf5/H5CX.c
		hdf5/H5D.c
		hdf5/H5Dbtree.c
		hdf5/H5Dbtree2.c
		hdf5/H5Dchunk.c
		hdf5/H5Dcompact.c
		hdf5/H5Dcontig.c
		hdf5/H5Ddbg.c
		hdf5/H5Ddeprec.c
		hdf5/H5Dearray.c
		hdf5/H5Defl.c
		hdf5/H5Df.c
		hdf5/H5Dfarray.c
		hdf5/H5Dfill.c
		hdf5/H5Dint.c
		hdf5/H5Dio.c
		hdf5/H5Dlayout.c
		hdf5/H5Dmpio.c
		hdf5/H5Dnone.c
		hdf5/H5Doh.c
		hdf5/H5Dscatgath.c
		hdf5/H5Dselect.c
		hdf5/H5Dsingle.c
		hdf5/H5Dtest.c
		hdf5/H5Dvirtual.c
		hdf5/H5E.c
		hdf5/H5EA.c
		hdf5/H5EAcache.c
		hdf5/H5EAdbg.c
		hdf5/H5EAdblkpage.c
		hdf5/H5EAdblock.c
		hdf5/H5EAhdr.c
		hdf5/H5EAiblock.c
		hdf5/H5EAint.c
		hdf5/H5EAsblock.c
		hdf5/H5EAstat.c
		hdf5/H5EAtest.c
		hdf5/H5Edeprec.c
		hdf5/H5Ef.c
		hdf5/H5Eint.c
		hdf5/H5F.c
		hdf5/H5FA.c
		hdf5/H5FAcache.c
		hdf5/H5FAdbg.c
		hdf5/H5FAdblkpage.c
		hdf5/H5FAdblock.c
		hdf5/H5FAhdr.c
		hdf5/H5FAint.c
		hdf5/H5FAstat.c
		hdf5/H5FAtest.c
		hdf5/H5FD.c
		hdf5/H5FDcore.c
		hdf5/H5FDdirect.c
		hdf5/H5FDfamily.c
		hdf5/H5FDint.c
		hdf5/H5FDlog.c
		hdf5/H5FDmpi.c
		hdf5/H5FDmpio.c
		hdf5/H5FDmulti.c
		hdf5/H5FDsec2.c
		hdf5/H5FDspace.c
		hdf5/H5FDstdio.c
		hdf5/H5FDtest.c
		hdf5/H5FDwindows.c
		hdf5/H5FL.c
		hdf5/H5FO.c
		hdf5/H5FS.c
		hdf5/H5FScache.c
		hdf5/H5FSdbg.c
		hdf5/H5FSint.c
		hdf5/H5FSsection.c
		hdf5/H5FSstat.c
		hdf5/H5FStest.c
		hdf5/H5Faccum.c
		hdf5/H5Fcwfs.c
		hdf5/H5Fdbg.c
		hdf5/H5Fdeprec.c
		hdf5/H5Fefc.c
		hdf5/H5Ff.c
		hdf5/H5Ffake.c
		hdf5/H5Fint.c
		hdf5/H5Fio.c
		hdf5/H5Fmount.c
		hdf5/H5Fmpi.c
		hdf5/H5Fquery.c
		hdf5/H5Fsfile.c
		hdf5/H5Fspace.c
		hdf5/H5Fsuper.c
		hdf5/H5Fsuper_cache.c
		hdf5/H5Ftest.c
		hdf5/H5G.c
		hdf5/H5Gbtree2.c
		hdf5/H5Gcache.c
		hdf5/H5Gcompact.c
		hdf5/H5Gdense.c
		hdf5/H5Gdeprec.c
		hdf5/H5Gent.c
		hdf5/H5Gf.c
		hdf5/H5Gint.c
		hdf5/H5Glink.c
		hdf5/H5Gloc.c
		hdf5/H5Gname.c
		hdf5/H5Gnode.c
		hdf5/H5Gobj.c
		hdf5/H5Goh.c
		hdf5/H5Groot.c
		hdf5/H5Gstab.c
		hdf5/H5Gtest.c
		hdf5/H5Gtraverse.c
		hdf5/H5HF.c
		hdf5/H5HFbtree2.c
		hdf5/H5HFcache.c
		hdf5/H5HFdbg.c
		hdf5/H5HFdblock.c
		hdf5/H5HFdtable.c
		hdf5/H5HFhdr.c
		hdf5/H5HFhuge.c
		hdf5/H5HFiblock.c
		hdf5/H5HFiter.c
		hdf5/H5HFman.c
		hdf5/H5HFsection.c
		hdf5/H5HFspace.c
		hdf5/H5HFstat.c
		hdf5/H5HFtest.c
		hdf5/H5HFtiny.c
		hdf5/H5HG.c
		hdf5/H5HGcache.c
		hdf5/H5HGdbg.c
		hdf5/H5HGquery.c
		hdf5/H5HL.c
		hdf5/H5HLcache.c
		hdf5/H5HLdbg.c
		hdf5/H5HLdblk.c
		hdf5/H5HLint.c
		hdf5/H5HLprfx.c
		hdf5/H5HP.c
		hdf5/H5I.c
		hdf5/H5If.c
		hdf5/H5Itest.c
		hdf5/H5L.c
		hdf5/H5Lexternal.c
		hdf5/H5lib_settings.c
		hdf5/H5Lf.c
		hdf5/H5MF.c
		hdf5/H5MFaggr.c
		hdf5/H5MFdbg.c
		hdf5/H5MFsection.c
		hdf5/H5MM.c
		hdf5/H5MP.c
		hdf5/H5MPtest.c
		hdf5/H5O.c
		hdf5/H5Oainfo.c
		hdf5/H5Oalloc.c
		hdf5/H5Oattr.c
		hdf5/H5Oattribute.c
		hdf5/H5Obogus.c
		hdf5/H5Obtreek.c
		hdf5/H5Ocache.c
		hdf5/H5Ocache_image.c
		hdf5/H5Ochunk.c
		hdf5/H5Ocont.c
		hdf5/H5Ocopy.c
		hdf5/H5Odbg.c
		hdf5/H5Odrvinfo.c
		hdf5/H5Odtype.c
		hdf5/H5Oefl.c
		hdf5/H5Of.c
		hdf5/H5Ofill.c
		hdf5/H5Oflush.c
		hdf5/H5Ofsinfo.c
		hdf5/H5Oginfo.c
		hdf5/H5Oint.c
		hdf5/H5Olayout.c
		hdf5/H5Olinfo.c
		hdf5/H5Olink.c
		hdf5/H5Omessage.c
		hdf5/H5Omtime.c
		hdf5/H5Oname.c
		hdf5/H5Onull.c
		hdf5/H5Opline.c
		hdf5/H5Orefcount.c
		hdf5/H5Osdspace.c
		hdf5/H5Oshared.c
		hdf5/H5Oshmesg.c
		hdf5/H5Ostab.c
		hdf5/H5Otest.c
		hdf5/H5Ounknown.c
		hdf5/H5P.c
		hdf5/H5PB.c
		hdf5/H5PL.c
		hdf5/H5PLint.c
		hdf5/H5PLpath.c
		hdf5/H5PLplugin_cache.c
		hdf5/H5Pacpl.c
		hdf5/H5Pdapl.c
		hdf5/H5Pdcpl.c
		hdf5/H5Pdeprec.c
		hdf5/H5Pdxpl.c
		hdf5/H5Pencdec.c
		hdf5/H5Pf.c
		hdf5/H5Pfapl.c
		hdf5/H5Pfcpl.c
		hdf5/H5Pfmpl.c
		hdf5/H5Pgcpl.c
		hdf5/H5Pint.c
		hdf5/H5Plapl.c
		hdf5/H5Plcpl.c
		hdf5/H5Pocpl.c
		hdf5/H5Pocpypl.c
		hdf5/H5Pstrcpl.c
		hdf5/H5Ptest.c
		hdf5/H5R.c
		hdf5/H5RS.c
		hdf5/H5Rdeprec.c
		hdf5/H5Rf.c
		hdf5/H5Rint.c
		hdf5/H5S.c
		hdf5/H5SL.c
		hdf5/H5SM.c
		hdf5/H5SMbtree2.c
		hdf5/H5SMcache.c
		hdf5/H5SMmessage.c
		hdf5/H5SMtest.c
		hdf5/H5ST.c
		hdf5/H5Sall.c
		hdf5/H5Sdbg.c
		hdf5/H5Sf.c
		hdf5/H5Shyper.c
		hdf5/H5Smpio.c
		hdf5/H5Snone.c
		hdf5/H5Spoint.c
		hdf5/H5Sselect.c
		hdf5/H5Stest.c
		hdf5/H5T.c
		hdf5/H5TS.c
		hdf5/H5Tarray.c
		hdf5/H5Tbit.c
		hdf5/H5Tcommit.c
		hdf5/H5Tcompound.c
		hdf5/H5Tconv.c
		hdf5/H5Tcset.c
		hdf5/H5Tdbg.c
		hdf5/H5Tdeprec.c
		hdf5/H5Tenum.c
		hdf5/H5Tf.c
		hdf5/H5Tfields.c
		hdf5/H5Tfixed.c
		hdf5/H5Tfloat.c
		hdf5/H5Tinit.c
		hdf5/H5Tnative.c
		hdf5/H5Toffset.c
		hdf5/H5Toh.c
		hdf5/H5Topaque.c
		hdf5/H5Torder.c
		hdf5/H5Tpad.c
		hdf5/H5Tprecis.c
		hdf5/H5Tstrpad.c
		hdf5/H5Tvisit.c
		hdf5/H5Tvlen.c
		hdf5/H5UC.c
		hdf5/H5VM.c
		hdf5/H5WB.c
		hdf5/H5Z.c
		hdf5/H5Zdeflate.c
		hdf5/H5Zf.c
		hdf5/H5Zfletcher32.c
		hdf5/H5Znbit.c
		hdf5/H5Zscaleoffset.c
		hdf5/H5Zshuffle.c
		hdf5/H5Ztrans.c
		hdf5/H5_f.c
		hdf5/H5checksum.c
		hdf5/H5dbg.c
		hdf5/H5f90kit.c
		hdf5/H5system.c
		hdf5/H5timer.c
		hdf5/H5trace.c
		hdf5/H5_gen.F90
		hdf5/H5Af.c
		hdf5/H5Aff.F90
		hdf5/H5Df.c
		hdf5/H5Dff.F90
		hdf5/H5Ef.c
		hdf5/H5Eff.F90
		hdf5/H5Ff.c
		hdf5/H5Fff.F90
		hdf5/H5Gf.c
		hdf5/H5Gff.F90
		hdf5/H5If.c
		hdf5/H5Iff.F90
		hdf5/H5Lf.c
		hdf5/H5Lff.F90
		hdf5/H5Of.c
		hdf5/H5Off.F90
		hdf5/H5Pf.c
		hdf5/H5Pff.F90
		hdf5/H5Rf.c
		hdf5/H5Rff.F90
		hdf5/H5Sf.c
		hdf5/H5Sff.F90
		hdf5/H5Tf.c
		hdf5/H5Tff.F90
		hdf5/H5Zf.c
		hdf5/H5Zff.F90
		hdf5/H5_f.c
		hdf5/H5_ff.F90
		hdf5/H5f90global.F90
		hdf5/H5f90kit.c
		hdf5/H5fortkit.F90
		hdf5/H5fortran_types.F90
		hdf5/HDF5.F90
		hdf5/H5DO.c
		hdf5/H5DS.c
		hdf5/H5IM.c
		hdf5/H5LD.c
		hdf5/H5LT.c
		hdf5/H5LTanalyze.c
		hdf5/H5LTparse.c
		hdf5/H5PT.c
		hdf5/H5TB.c
		hdf5/H5DSfc.c
		hdf5/H5DSff.F90
		hdf5/H5IMcc.c
		hdf5/H5IMfc.c
		hdf5/H5IMff.F90
		hdf5/H5LTfc.c
		hdf5/H5LTff.F90
		hdf5/H5TBfc.c
		hdf5/H5TBff.F90
		libz/adler32.c
		libz/crc32.c
		libz/deflate.c
		libz/infback.c
		libz/inffast.c
		libz/inflate.c
		libz/inftrees.c
		libz/trees.c
		libz/zutil.c
		libz/compress.c
		libz/uncompr.c
		libz/gzclose.c
		libz/gzlib.c
		libz/gzread.c
		libz/gzwrite.c
	     );

$need_math = 'no';

$bench_flags = '';
$bench_fflags = '';

sub invoke {
	my ($me) = @_;
	my @rc;
        my @temp = main::read_file('control');
	my $exe = $me->exe_file;
        for my $line (@temp) {
            my ($runid) = split(/\s+/, $line);
	    push (@rc, { 'command' => $exe, 
			'args'    => [$runid], 
			'output'  => "$exename.stdout.out",
			'error'   => "$exename.stderr.err",
			});
        }
	return @rc;
}

%deps = (
		'pot3d.F90' => [
		'number_types.F90',
		'zm_parse_modules.F90',
		'zm_parse.F90',
		'zm_sds_modules.F90',
		'zm_sds.F90',
		'hdf5/HDF5.F90',
		'specmpitime/specmpitime_mod.F90'
		],
		'zm_parse_modules.F90' => [
		'number_types.F90',
		'hdf5/HDF5.F90'
		],
		'zm_parse.F90' => [
		'number_types.F90',
		'zm_parse_modules.F90',
		'hdf5/HDF5.F90'
		],
		'zm_sds_modules.F90' => [
		'number_types.F90',
		'hdf5/HDF5.F90'
		],
		'zm_sds.F90' => [
		'number_types.F90',
		'zm_sds_modules.F90',
		'hdf5/HDF5.F90',
		'hdf5/H5DSff.F90'
		],
	   'hdf5/HDF5.F90' => [
		'hdf5/H5f90global.F90',
		'hdf5/H5_ff.F90',
		'hdf5/H5_gen.F90',
		'hdf5/H5Zff.F90',
		'hdf5/H5Rff.F90',
		'hdf5/H5Pff.F90',
		'hdf5/H5Off.F90',
		'hdf5/H5Tff.F90',
		'hdf5/H5Aff.F90',
		'hdf5/H5Dff.F90',
		'hdf5/H5Sff.F90',
		'hdf5/H5Lff.F90',
		'hdf5/H5Iff.F90',
		'hdf5/H5Eff.F90',
		'hdf5/H5Gff.F90',
		'hdf5/H5Fff.F90'
		],
	   'hdf5/H5f90global.F90' => [
		'hdf5/H5fortran_types.F90'
		],
	   'hdf5/H5_gen.F90' => [
		'hdf5/H5f90global.F90',
		'hdf5/H5Aff.F90',
		'hdf5/H5Dff.F90',
		'hdf5/H5Pff.F90'
		],
	   'hdf5/H5_ff.F90' => [
		'hdf5/H5f90global.F90',
		'hdf5/H5Fff.F90',
		],
	   'hdf5/H5Aff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Dff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Eff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Fff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Gff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Iff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Lff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Off.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Pff.F90' => [
		'hdf5/H5f90global.F90',
		'hdf5/H5fortkit.F90'
		],
	   'hdf5/H5Rff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Sff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Tff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5Zff.F90' => [
		'hdf5/H5f90global.F90'
		],
	   'hdf5/H5DSff.F90' => [
		'hdf5/H5fortran_types.F90',
		'hdf5/HDF5.F90'
		],
	   'hdf5/H5IMff.F90' => [
		'hdf5/H5fortran_types.F90',
		'hdf5/HDF5.F90'
		],
	   'hdf5/H5LTff.F90' => [
		'hdf5/H5fortran_types.F90',
		'hdf5/HDF5.F90'
		],
	   'hdf5/H5TBff.F90' => [
		'hdf5/H5fortran_types.F90',
		'hdf5/HDF5.F90'
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
