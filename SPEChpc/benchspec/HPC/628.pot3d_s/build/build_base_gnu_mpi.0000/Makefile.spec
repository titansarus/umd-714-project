TUNE=base
LABEL=gnu_mpi
NUMBER=628
NAME=pot3d_s
SOURCES= number_types.F90 pot3d.F90 zm_parse.F90 zm_parse_modules.F90 \
	 zm_sds.F90 zm_sds_modules.F90 specmpitime/specmpitime_mod.F90 \
	 specmpitime/specmpitime.c hdf5/H5.c hdf5/H5A.c hdf5/H5AC.c \
	 hdf5/H5ACdbg.c hdf5/H5ACmpio.c hdf5/H5ACproxy_entry.c hdf5/H5Abtree2.c \
	 hdf5/H5Adense.c hdf5/H5Adeprec.c hdf5/H5Af.c hdf5/H5Aint.c \
	 hdf5/H5Atest.c hdf5/H5B.c hdf5/H5B2.c hdf5/H5B2cache.c hdf5/H5B2dbg.c \
	 hdf5/H5B2hdr.c hdf5/H5B2int.c hdf5/H5B2internal.c hdf5/H5B2leaf.c \
	 hdf5/H5B2stat.c hdf5/H5B2test.c hdf5/H5Bcache.c hdf5/H5Bdbg.c hdf5/H5C.c \
	 hdf5/H5CS.c hdf5/H5Cdbg.c hdf5/H5Cepoch.c hdf5/H5Cimage.c hdf5/H5Clog.c \
	 hdf5/H5Clog_json.c hdf5/H5Clog_trace.c hdf5/H5Cmpio.c \
	 hdf5/H5Cprefetched.c hdf5/H5Cquery.c hdf5/H5Ctag.c hdf5/H5Ctest.c \
	 hdf5/H5CX.c hdf5/H5D.c hdf5/H5Dbtree.c hdf5/H5Dbtree2.c hdf5/H5Dchunk.c \
	 hdf5/H5Dcompact.c hdf5/H5Dcontig.c hdf5/H5Ddbg.c hdf5/H5Ddeprec.c \
	 hdf5/H5Dearray.c hdf5/H5Defl.c hdf5/H5Df.c hdf5/H5Dfarray.c \
	 hdf5/H5Dfill.c hdf5/H5Dint.c hdf5/H5Dio.c hdf5/H5Dlayout.c \
	 hdf5/H5Dmpio.c hdf5/H5Dnone.c hdf5/H5Doh.c hdf5/H5Dscatgath.c \
	 hdf5/H5Dselect.c hdf5/H5Dsingle.c hdf5/H5Dtest.c hdf5/H5Dvirtual.c \
	 hdf5/H5E.c hdf5/H5EA.c hdf5/H5EAcache.c hdf5/H5EAdbg.c \
	 hdf5/H5EAdblkpage.c hdf5/H5EAdblock.c hdf5/H5EAhdr.c hdf5/H5EAiblock.c \
	 hdf5/H5EAint.c hdf5/H5EAsblock.c hdf5/H5EAstat.c hdf5/H5EAtest.c \
	 hdf5/H5Edeprec.c hdf5/H5Ef.c hdf5/H5Eint.c hdf5/H5F.c hdf5/H5FA.c \
	 hdf5/H5FAcache.c hdf5/H5FAdbg.c hdf5/H5FAdblkpage.c hdf5/H5FAdblock.c \
	 hdf5/H5FAhdr.c hdf5/H5FAint.c hdf5/H5FAstat.c hdf5/H5FAtest.c \
	 hdf5/H5FD.c hdf5/H5FDcore.c hdf5/H5FDdirect.c hdf5/H5FDfamily.c \
	 hdf5/H5FDint.c hdf5/H5FDlog.c hdf5/H5FDmpi.c hdf5/H5FDmpio.c \
	 hdf5/H5FDmulti.c hdf5/H5FDsec2.c hdf5/H5FDspace.c hdf5/H5FDstdio.c \
	 hdf5/H5FDtest.c hdf5/H5FDwindows.c hdf5/H5FL.c hdf5/H5FO.c hdf5/H5FS.c \
	 hdf5/H5FScache.c hdf5/H5FSdbg.c hdf5/H5FSint.c hdf5/H5FSsection.c \
	 hdf5/H5FSstat.c hdf5/H5FStest.c hdf5/H5Faccum.c hdf5/H5Fcwfs.c \
	 hdf5/H5Fdbg.c hdf5/H5Fdeprec.c hdf5/H5Fefc.c hdf5/H5Ff.c hdf5/H5Ffake.c \
	 hdf5/H5Fint.c hdf5/H5Fio.c hdf5/H5Fmount.c hdf5/H5Fmpi.c hdf5/H5Fquery.c \
	 hdf5/H5Fsfile.c hdf5/H5Fspace.c hdf5/H5Fsuper.c hdf5/H5Fsuper_cache.c \
	 hdf5/H5Ftest.c hdf5/H5G.c hdf5/H5Gbtree2.c hdf5/H5Gcache.c \
	 hdf5/H5Gcompact.c hdf5/H5Gdense.c hdf5/H5Gdeprec.c hdf5/H5Gent.c \
	 hdf5/H5Gf.c hdf5/H5Gint.c hdf5/H5Glink.c hdf5/H5Gloc.c hdf5/H5Gname.c \
	 hdf5/H5Gnode.c hdf5/H5Gobj.c hdf5/H5Goh.c hdf5/H5Groot.c hdf5/H5Gstab.c \
	 hdf5/H5Gtest.c hdf5/H5Gtraverse.c hdf5/H5HF.c hdf5/H5HFbtree2.c \
	 hdf5/H5HFcache.c hdf5/H5HFdbg.c hdf5/H5HFdblock.c hdf5/H5HFdtable.c \
	 hdf5/H5HFhdr.c hdf5/H5HFhuge.c hdf5/H5HFiblock.c hdf5/H5HFiter.c \
	 hdf5/H5HFman.c hdf5/H5HFsection.c hdf5/H5HFspace.c hdf5/H5HFstat.c \
	 hdf5/H5HFtest.c hdf5/H5HFtiny.c hdf5/H5HG.c hdf5/H5HGcache.c \
	 hdf5/H5HGdbg.c hdf5/H5HGquery.c hdf5/H5HL.c hdf5/H5HLcache.c \
	 hdf5/H5HLdbg.c hdf5/H5HLdblk.c hdf5/H5HLint.c hdf5/H5HLprfx.c \
	 hdf5/H5HP.c hdf5/H5I.c hdf5/H5If.c hdf5/H5Itest.c hdf5/H5L.c \
	 hdf5/H5Lexternal.c hdf5/H5lib_settings.c hdf5/H5Lf.c hdf5/H5MF.c \
	 hdf5/H5MFaggr.c hdf5/H5MFdbg.c hdf5/H5MFsection.c hdf5/H5MM.c \
	 hdf5/H5MP.c hdf5/H5MPtest.c hdf5/H5O.c hdf5/H5Oainfo.c hdf5/H5Oalloc.c \
	 hdf5/H5Oattr.c hdf5/H5Oattribute.c hdf5/H5Obogus.c hdf5/H5Obtreek.c \
	 hdf5/H5Ocache.c hdf5/H5Ocache_image.c hdf5/H5Ochunk.c hdf5/H5Ocont.c \
	 hdf5/H5Ocopy.c hdf5/H5Odbg.c hdf5/H5Odrvinfo.c hdf5/H5Odtype.c \
	 hdf5/H5Oefl.c hdf5/H5Of.c hdf5/H5Ofill.c hdf5/H5Oflush.c \
	 hdf5/H5Ofsinfo.c hdf5/H5Oginfo.c hdf5/H5Oint.c hdf5/H5Olayout.c \
	 hdf5/H5Olinfo.c hdf5/H5Olink.c hdf5/H5Omessage.c hdf5/H5Omtime.c \
	 hdf5/H5Oname.c hdf5/H5Onull.c hdf5/H5Opline.c hdf5/H5Orefcount.c \
	 hdf5/H5Osdspace.c hdf5/H5Oshared.c hdf5/H5Oshmesg.c hdf5/H5Ostab.c \
	 hdf5/H5Otest.c hdf5/H5Ounknown.c hdf5/H5P.c hdf5/H5PB.c hdf5/H5PL.c \
	 hdf5/H5PLint.c hdf5/H5PLpath.c hdf5/H5PLplugin_cache.c hdf5/H5Pacpl.c \
	 hdf5/H5Pdapl.c hdf5/H5Pdcpl.c hdf5/H5Pdeprec.c hdf5/H5Pdxpl.c \
	 hdf5/H5Pencdec.c hdf5/H5Pf.c hdf5/H5Pfapl.c hdf5/H5Pfcpl.c \
	 hdf5/H5Pfmpl.c hdf5/H5Pgcpl.c hdf5/H5Pint.c hdf5/H5Plapl.c \
	 hdf5/H5Plcpl.c hdf5/H5Pocpl.c hdf5/H5Pocpypl.c hdf5/H5Pstrcpl.c \
	 hdf5/H5Ptest.c hdf5/H5R.c hdf5/H5RS.c hdf5/H5Rdeprec.c hdf5/H5Rf.c \
	 hdf5/H5Rint.c hdf5/H5S.c hdf5/H5SL.c hdf5/H5SM.c hdf5/H5SMbtree2.c \
	 hdf5/H5SMcache.c hdf5/H5SMmessage.c hdf5/H5SMtest.c hdf5/H5ST.c \
	 hdf5/H5Sall.c hdf5/H5Sdbg.c hdf5/H5Sf.c hdf5/H5Shyper.c hdf5/H5Smpio.c \
	 hdf5/H5Snone.c hdf5/H5Spoint.c hdf5/H5Sselect.c hdf5/H5Stest.c \
	 hdf5/H5T.c hdf5/H5TS.c hdf5/H5Tarray.c hdf5/H5Tbit.c hdf5/H5Tcommit.c \
	 hdf5/H5Tcompound.c hdf5/H5Tconv.c hdf5/H5Tcset.c hdf5/H5Tdbg.c \
	 hdf5/H5Tdeprec.c hdf5/H5Tenum.c hdf5/H5Tf.c hdf5/H5Tfields.c \
	 hdf5/H5Tfixed.c hdf5/H5Tfloat.c hdf5/H5Tinit.c hdf5/H5Tnative.c \
	 hdf5/H5Toffset.c hdf5/H5Toh.c hdf5/H5Topaque.c hdf5/H5Torder.c \
	 hdf5/H5Tpad.c hdf5/H5Tprecis.c hdf5/H5Tstrpad.c hdf5/H5Tvisit.c \
	 hdf5/H5Tvlen.c hdf5/H5UC.c hdf5/H5VM.c hdf5/H5WB.c hdf5/H5Z.c \
	 hdf5/H5Zdeflate.c hdf5/H5Zf.c hdf5/H5Zfletcher32.c hdf5/H5Znbit.c \
	 hdf5/H5Zscaleoffset.c hdf5/H5Zshuffle.c hdf5/H5Ztrans.c hdf5/H5_f.c \
	 hdf5/H5checksum.c hdf5/H5dbg.c hdf5/H5f90kit.c hdf5/H5system.c \
	 hdf5/H5timer.c hdf5/H5trace.c hdf5/H5_gen.F90 hdf5/H5Af.c hdf5/H5Aff.F90 \
	 hdf5/H5Df.c hdf5/H5Dff.F90 hdf5/H5Ef.c hdf5/H5Eff.F90 hdf5/H5Ff.c \
	 hdf5/H5Fff.F90 hdf5/H5Gf.c hdf5/H5Gff.F90 hdf5/H5If.c hdf5/H5Iff.F90 \
	 hdf5/H5Lf.c hdf5/H5Lff.F90 hdf5/H5Of.c hdf5/H5Off.F90 hdf5/H5Pf.c \
	 hdf5/H5Pff.F90 hdf5/H5Rf.c hdf5/H5Rff.F90 hdf5/H5Sf.c hdf5/H5Sff.F90 \
	 hdf5/H5Tf.c hdf5/H5Tff.F90 hdf5/H5Zf.c hdf5/H5Zff.F90 hdf5/H5_f.c \
	 hdf5/H5_ff.F90 hdf5/H5f90global.F90 hdf5/H5f90kit.c hdf5/H5fortkit.F90 \
	 hdf5/H5fortran_types.F90 hdf5/HDF5.F90 hdf5/H5DO.c hdf5/H5DS.c \
	 hdf5/H5IM.c hdf5/H5LD.c hdf5/H5LT.c hdf5/H5LTanalyze.c hdf5/H5LTparse.c \
	 hdf5/H5PT.c hdf5/H5TB.c hdf5/H5DSfc.c hdf5/H5DSff.F90 hdf5/H5IMcc.c \
	 hdf5/H5IMfc.c hdf5/H5IMff.F90 hdf5/H5LTfc.c hdf5/H5LTff.F90 \
	 hdf5/H5TBfc.c hdf5/H5TBff.F90 libz/adler32.c libz/crc32.c libz/deflate.c \
	 libz/infback.c libz/inffast.c libz/inflate.c libz/inftrees.c \
	 libz/trees.c libz/zutil.c libz/compress.c libz/uncompr.c libz/gzclose.c \
	 libz/gzlib.c libz/gzread.c libz/gzwrite.c
EXEBASE=pot3d
NEED_MATH=
USE_SPACK=
BENCHLANG=F

BENCH_CFLAGS     = -Ihdf5/include -Ilibz/include -DNDEBUG -Ispecmpitime 
BENCH_FFLAGS     = 
BENCH_FLAGS      = 
BENCH_FPPFLAGS   = -Ihdf5  -Ihdf5/include -DNDEBUG 
CC               = mpiicc
CC_VERSION_OPTION = --version
COPTIMIZE        = -Ofast -march=native -lm
CXX              = mpiicpc
CXXOPTIMIZE      = -Ofast -march=native -std=c++14
CXX_VERSION_OPTION = --version
FC               = mpiifort
FC_VERSION_OPTION = --version
FOPTIMIZE        = -Ofast -march=native -fno-stack-protector
FPORTABILITY     = -ffree-line-length-none
OS               = unix
PMODEL           = 
absolutely_no_locking = 0
abstol           = 4e-07
action           = build
allow_label_override = 0
backup_config    = 1
baseexe          = pot3d
basepeak         = 0
benchargs        = 
benchdir         = benchspec
benchmark        = 628.pot3d_s
binary           = 
bindir           = exe
builddir         = build
bundleaction     = 
bundlename       = 
calctol          = 1
changedhash      = 0
check_version    = 0
clean_between_builds = no
command_add_redirect = 0
commanderrfile   = speccmds.err
commandexe       = pot3d_base.gnu_mpi
commandfile      = speccmds.cmd
commandoutfile   = speccmds.out
commandstdoutfile = speccmds.stdout
comparedir       = compare
compareerrfile   = compare.err
comparefile      = compare.cmd
compareoutfile   = compare.out
comparestdoutfile = compare.stdout
compile_error    = 0
compwhite        = 
configdir        = config
configfile       = Project-try1
configpath       = /home/jumeike/scratch.cmsc714/final-project/config/Project-try1.cfg
copies           = 1
current_range    = 
datadir          = data
default_size     = ref
default_submit   = $command
delay            = 0
deletebinaries   = 0
deletework       = 0
dependent_workloads = 0
device           = 
difflines        = 10
dirprot          = 511
discard_power_samples = 1
enable_monitor   = 1
endian           = 12345678
env_vars         = 0
expand_notes     = 0
expid            = 
exthash_bits     = 0
failflags        = 0
fake             = 0
feedback         = 1
flag_url_base    = https://www.spec.org/auto/hpc2021/Docs/benchmarks/flags/
floatcompare     = 
force_monitor    = 0
hostname         = login-1.zaratan.umd.edu
http_proxy       = 
http_timeout     = 30
hw_avail         = Nov-2099
idle_current_range = 
idledelay        = 10
idleduration     = 60
ignore_errors    = 0
ignore_sigint    = 0
ignorecase       = 
info_wrap_columns = 50
inputdir         = input
inputgenerrfile  = inputgen.err
inputgenfile     = inputgen.cmd
inputgenoutfile  = inputgen.out
inputgenstdoutfile = inputgen.stdout
interconnect_fs_hw_model = BI 100 Series
interconnect_fs_hw_switch_fs_count = 1
interconnect_fs_hw_switch_fs_data_rate = 100 Gb/s
interconnect_fs_hw_switch_fs_firmware = 10.3.0.0.60
interconnect_fs_hw_switch_fs_model000 = BI 100 Series 48 Port 2
interconnect_fs_hw_switch_fs_model001 = PSU
interconnect_fs_hw_switch_fs_ports = 48
interconnect_fs_hw_topo = Mesh
interconnect_fs_hw_vendor = Big Interconnect Company
interconnect_fs_order = 0
interconnect_fs_purpose = MPI Traffic
interconnect_fs_syslbl = Big Interconnect Company
islibrary        = 0
iteration        = -1
iterations       = 2
keeptmp          = 0
label            = gnu_mpi
license_num      = 9999
line_width       = 0
link_input_files = 1
locking          = 1
log              = hpc2021
log_line_width   = 0
log_timestamp    = 0
logname          = /home/jumeike/scratch.cmsc714/final-project/result/hpc2021.064.log
lognum           = 064
mail_reports     = all
mailcompress     = 0
mailmethod       = smtp
mailport         = 25
mailserver       = 127.0.0.1
mailto           = 
make             = specmake
make_no_clobber  = 0
makefile_template = Makefile.YYYtArGeTYYYspec
makeflags        = -j 40
max_average_uncertainty = 1
max_hum_limit    = 0
max_report_runs  = 0
max_unknown_uncertainty = 1
mean_anyway      = 0
meter_connect_timeout = 30
meter_errors_default = 5
meter_errors_percentage = 5
min_report_runs  = 2
min_temp_limit   = 20
minimize_builddirs = 0
minimize_rundirs = 0
name             = pot3d_s
nansupport       = 
need_math        = 
no_input_handler = close
no_monitor       = 
node_compute_count = 2
node_compute_hw_accel_connect = PCIe 3.0 16x
node_compute_hw_accel_count = 4
node_compute_hw_accel_desc = See Notes
node_compute_hw_accel_ecc = Yes
node_compute_hw_accel_model = Tesla V100-PCIE-16GB
node_compute_hw_accel_type = GPU
node_compute_hw_accel_vendor = NVIDIA Corporation
node_compute_hw_adapter_fs_count = 0
node_compute_hw_adapter_fs_data_rate = None
node_compute_hw_adapter_fs_driver = None
node_compute_hw_adapter_fs_firmware = None
node_compute_hw_adapter_fs_interconnect = None
node_compute_hw_adapter_fs_model = None
node_compute_hw_adapter_fs_ports_used = 0
node_compute_hw_adapter_fs_slot_type = None
node_compute_hw_cpu_char = Turbo up to 3.4 GHz
node_compute_hw_cpu_mhz = 2250
node_compute_hw_cpu_name = Turbo CPU
node_compute_hw_disk = 1 x 480 GB  SATA 2.5" SSD
node_compute_hw_memory = 256 GB (8 x 32 GB 2Rx8 PC4-3200AA-R)
node_compute_hw_model = Turblaster 5000
node_compute_hw_nchips = 1
node_compute_hw_ncores = 128
node_compute_hw_ncoresperchip = 64
node_compute_hw_ncpuorder = 1 chips
node_compute_hw_nthreadspercore = 1
node_compute_hw_ocache = None
node_compute_hw_other = None
node_compute_hw_pcache = 32 KB I + 32 KB D on chip per core
node_compute_hw_scache = 512 KB I+D on chip per core
node_compute_hw_tcache000 = 256 MB I+D on chip per chip
node_compute_hw_tcache001 = 16 MB shared / 4 cores
node_compute_hw_vendor = Mega Technology
node_compute_order = 1
node_compute_purpose = compute
node_compute_sw_localfile = xfs
node_compute_sw_os000 = SUSE Linux Enterprise Linux Server 12
node_compute_sw_os001 = 4.12.14-94.41-default
node_compute_sw_other = None
node_compute_sw_sharedfile = None
node_compute_sw_state = Multi-user, run level 3
node_compute_syslbl = TurboBlaster 5000
node_fileserver_count = 1
node_fileserver_hw_adapter_fs_count = 1
node_fileserver_hw_adapter_fs_data_rate = 100 Gb/s
node_fileserver_hw_adapter_fs_driver = 10.9.1.0.15
node_fileserver_hw_adapter_fs_firmware = 10.9.0.1.0
node_fileserver_hw_adapter_fs_interconnect = BG 5000 series
node_fileserver_hw_adapter_fs_ports_used = 1
node_fileserver_hw_adapter_fs_slot_type = PCI-Express 3.0 x16
node_fileserver_hw_cpu_char = None
node_fileserver_hw_cpu_mhz = 2700
node_fileserver_hw_cpu_name = Intel Xeon Platinum 8280
node_fileserver_hw_disk = 1 x 1 TB 12 Gbps SAS 2.5" SSD (JBOD)
node_fileserver_hw_memory = 768 GB (24 x 32 GB 2Rx4 PC4-2933Y-R)
node_fileserver_hw_model = BG650
node_fileserver_hw_nchips = 2
node_fileserver_hw_ncores = 56
node_fileserver_hw_ncoresperchip = 28
node_fileserver_hw_ncpuorder = 1-2 chips
node_fileserver_hw_nthreadspercore = 1
node_fileserver_hw_ocache = None
node_fileserver_hw_other = None
node_fileserver_hw_pcache = 32 KB I + 32 KB D on chip per core
node_fileserver_hw_scache = 1 MB I+D on chip per core
node_fileserver_hw_tcache = 39424 KB I+D on chip per chip
node_fileserver_hw_vendor = Big Storage Company
node_fileserver_order = 2
node_fileserver_purpose = Fileserver
node_fileserver_sw_localfile = None
node_fileserver_sw_os = Red Hat Enterprise Linux Server release 7.6
node_fileserver_sw_other = None
node_fileserver_sw_sharedfile = NFS
node_fileserver_sw_state = Multi-User, run level 3
node_fileserver_syslbl = NFS
noratios         = 0
note_preenv      = 0
notes_000        =  Environment Settings:
notes_005        =   Any extra settings
notes_submit_000 =     mpirun -np $ranks $command
notes_wrap_columns = 0
notes_wrap_indent =   
num              = 628
obiwan           = 
oldhash          = 
os_exe_ext       = 
output_format    = text
output_root      = 
outputdir        = output
parallel_test    = 0
parallel_test_submit = 0
parallel_test_workloads = 
path             = /home/jumeike/scratch.cmsc714/final-project/benchspec/HPC/628.pot3d_s
plain_train      = 0
platform         = 
pmodel           = MPI
power            = 0
preenv           = 1
prefix           = 
ranks            = 1
rawhash_bits     = 256
rebuild          = 0
reftime          = reftime
reltol           = 0.006
reportable       = 0
resultdir        = result
review           = 0
run              = all
rundir           = run
runhpc           = /home/jumeike/scratch.cmsc714/final-project/bin/harness/runhpc --config=Project-try1 --action=build small
runmode          = speed
safe_eval        = 1
save_build_files = 
section_specifier_fatal = 1
setprocgroup     = 1
setup_error      = 0
showtimer        = 0
sigint           = 2
size             = ref
size_class       = ref
skipabstol       = 
skipobiwan       = 
skipreltol       = 
skiptol          = 
smarttune        = base
specdiff         = specdiff
specrun          = specinvoke
srcalt           = 
srcdir           = src
srcsource        = /home/jumeike/scratch.cmsc714/final-project/benchspec/HPC/628.pot3d_s/src
stagger          = 10
strict_rundir_verify = 1
submit_default   = mpirun ${MPIRUN_OPTS} -np $ranks $command
sw_avail         = Nov-2099
sw_compiler000   = C/C++/Fortran: Version 10.2 of
sw_compiler001   = GNU Compilers
sw_mpi_library   = OpenMPI Version 3.1.5
sw_mpi_other     = None
sw_other         = None
sysinfo_hash_bits = 256
sysinfo_program  = 
sysinfo_program_hash = 
system_name      = Big Compute
system_vendor    = Mega Technology
table            = 1
teeout           = yes
test_date        = May-2025
test_sponsor     = Sponsor Name
tester           = Testing Company Name
threads          = 1
timerfile        = spectimes.txt
top              = /home/jumeike/scratch.cmsc714/final-project
train_single_thread = 0
train_with       = train
tune             = base
uid              = 449016
unbuffer         = 1
uncertainty_exception = 5
update           = 0
update_url       = http://www.spec.org/auto/hpc2021/updates/
use_spack        = 
use_submit_for_compare = 0
use_submit_for_speed = 1
username         = jumeike
verbose          = 5
verify_binaries  = 1
version          = 1.001007
version_url      = http://www.spec.org/auto/hpc2021/current_version
voltage_range    = 
worklist         = list
OUTPUT_RMFILES   = pot3d1.out
