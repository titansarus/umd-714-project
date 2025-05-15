#pragma once

struct Phase;

//! \brief setup the OpenMP Target devices according to the commandline arguments.
//! \return Errorcode
int set_openmp_devices(const struct Phase*const p);

//! \brief setup the OpenACC devices according to the commandline arguments.
//! \param p System for which the devices are set.
//! \return Errorcode
int set_openacc_devices(const struct Phase*const p);

//! get selected target device
int soma_get_device();
