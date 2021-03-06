!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eig66_hdf
#include "juDFT_env.h"
   !*****************************************************************
   ! DESC:Module for hdf-io of eig-file
   !      To be compatible with f90 interface of HDF, use kind for vars
   !
   !      !ATTENTION before calling openeig and after calling closeeig!
   !      !the hdf library has to be initialized or finalized, respectively
   !
   !      CONTAINS the following subroutines:
   !      openeig        opens file
   !      closeeig       closes file
   !      read_keb       reads kpt, enpara and basis data
   !      read_neig      read no of eigenvalues (and eigenvalues itself)
   !      read_eig       reads eigenvectors
   !      writeeig       saves all data for kpt
   !      writesingleeig saves data for one kpt and energy
   !
   !
   !                          Daniel Wortmann
   !*****************************************************************
   USE m_eig66_data
   USE m_types
#ifdef CPP_HDF
   USE hdf5
   USE m_hdf_tools
   IMPLICIT NONE

   PRIVATE
   INTEGER, PARAMETER :: one = 1, two = 2, three = 3, zero = 0
   !to have the correct
   !type for array constructors

#endif
   PUBLIC open_eig, close_eig
   PUBLIC read_eig
   PUBLIC write_eig!,writesingleeig,writeeigc,writebas

CONTAINS
   SUBROUTINE priv_find_data(id, d)
      INTEGER, INTENT(IN)::id
      TYPE(t_data_hdf), POINTER:: d

      CLASS(t_data), POINTER   ::dp
      CALL eig66_find_data(dp, id)
      SELECT TYPE (dp)
      TYPE is (t_data_hdf)
         d => dp
      CLASS default
         CALL judft_error("BUG: wrong datatype in eig66_hdf")
      END SELECT
   END SUBROUTINE priv_find_data
   !----------------------------------------------------------------------
   SUBROUTINE open_eig(id, mpi_comm, nmat, neig, nkpts, jspins, create, l_real, l_soc, readonly, l_olap, filename)

      !*****************************************************************
      !     opens hdf-file for eigenvectors+values
      !*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: id, mpi_comm
      INTEGER, INTENT(IN) :: nmat, neig, nkpts, jspins
      LOGICAL, INTENT(IN) :: create, readonly, l_real, l_soc, l_olap
      CHARACTER(LEN=*), OPTIONAL :: filename

#ifdef CPP_HDF

      INTEGER         :: hdferr, access_mode
      INTEGER(HID_T)  :: creation_prp, access_prp, spaceid
      LOGICAL         :: l_exist
      INTEGER(HSIZE_T):: dims(7)
      TYPE(t_data_HDF), POINTER::d
      !Set creation and access properties

#ifdef CPP_HDFMPI
      INCLUDE 'mpif.h'
      IF (readonly) THEN
         access_prp = H5P_DEFAULT_f
         creation_prp = H5P_DEFAULT_f
      ELSE
         CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr)
         !      CALL h5pset_fapl_mpiposix_f(access_prp,MPI_COMM,
         !     +.false.,hdferr)
         CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM, MPI_INFO_NULL, hdferr)
         creation_prp = H5P_DEFAULT_f !no special creation property
      ENDIF
#else
      access_prp = H5P_DEFAULT_f
      creation_prp = H5P_DEFAULT_f
#endif
      
      if(l_olap) call juDFT_error("olap not implemented for hdf5")
      
      CALL priv_find_data(id, d)
      IF (PRESENT(filename)) d%fname = filename
      CALL eig66_data_storedefault(d, jspins, nkpts, nmat, neig, l_real, l_soc)
      !set access_flags according
      IF (readonly) THEN
         access_mode = H5F_ACC_RDONLY_F
      ELSE
         access_mode = H5F_ACC_RDWR_F
      ENDIF
      !     OPEN FILE and get D%FID's
      IF (create) THEN
         INQUIRE (FILE=TRIM(d%fname)//'.hdf', EXIST=l_exist)
         access_mode = H5F_ACC_TRUNC_F
         !         IF (l_exist) WRITE (*,*)'Warning: eig.hdf was overwritten'
         CALL h5fcreate_f(TRIM(d%fname)//'.hdf', access_Mode, d%fid, hdferr, creation_prp, access_prp)
         ! create dataspaces and datasets
         !   scalars
         dims(:2) = (/nkpts, jspins/)
         CALL h5screate_simple_f(2, dims(:2), spaceid, hdferr)
         CALL h5dcreate_f(d%fid, "neig", H5T_NATIVE_INTEGER, spaceid, d%neigsetid, hdferr)
         CALL h5sclose_f(spaceid, hdferr)
         !     ew
         dims(:3) = (/neig, nkpts, jspins/)
         CALL h5screate_simple_f(3, dims(:3), spaceid, hdferr)
         CALL h5dcreate_f(d%fid, "energy", H5T_NATIVE_DOUBLE, spaceid, d%energysetid, hdferr)
         CALL h5sclose_f(spaceid, hdferr)
         !     ev
         if (l_real .and. .not. l_soc) THEN
            dims(:5) = (/one, nmat, neig, nkpts, jspins/)
         else
            dims(:5) = (/two, nmat, neig, nkpts, jspins/)
         endif
         CALL h5screate_simple_f(5, dims(:5), spaceid, hdferr)
         CALL h5dcreate_f(d%fid, "ev", H5T_NATIVE_DOUBLE, spaceid, d%evsetid, hdferr)
         CALL h5sclose_f(spaceid, hdferr)
      ELSE
         CALL h5fopen_f(TRIM(d%fname)//'.hdf', access_Mode, d%fid, hdferr, access_prp)
         !get dataset-ids
         CALL h5dopen_f(d%fid, 'energy', d%energysetid, hdferr)
         CALL h5dopen_f(d%fid, 'neig', d%neigsetid, hdferr)
         CALL h5dopen_f(d%fid, 'ev', d%evsetid, hdferr)
      endif
      IF (.NOT. access_prp == H5P_DEFAULT_f) CALL H5Pclose_f(access_prp&
              &     , hdferr)
#else
      CALL juDFT_error("Could not use HDF5 for IO, please recompile")
#endif
   END SUBROUTINE open_eig
   !----------------------------------------------------------------------
   SUBROUTINE close_eig(id, filename)
      !*****************************************************************
      !     closes hdf-file for eigenvectors+values
      !*****************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN)                   :: id
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filename

      INTEGER::hdferr
      TYPE(t_data_HDF), POINTER::d

      !close datasets
#ifdef CPP_HDF
      CALL priv_find_data(id, d)

      CALL h5dclose_f(d%energysetid, hdferr)
      CALL h5dclose_f(d%wikssetid, hdferr)
      CALL h5dclose_f(d%neigsetid, hdferr)
      CALL h5dclose_f(d%evsetid, hdferr)
      !close file
      CALL h5fclose_f(d%fid, hdferr)
      !If a filename was given and the name is not the current filename
      IF (PRESENT(filename)) THEN
         IF (filename .NE. d%fname) THEN
            CALL system("mv "//TRIM(d%fname)//".hdf "//TRIM(filename)//".hdf")
         ENDIF
      ENDIF
      d%fname = "eig"
      CALL eig66_remove_data(id)

#endif
   END SUBROUTINE close_eig
#ifdef CPP_HDF
   !----------------------------------------------------------------------
   SUBROUTINE priv_r_vec(d, nk, jspin, list, z)

      USE m_hdf_tools
      IMPLICIT NONE
      TYPE(t_data_HDF), INTENT(IN)::d
      INTEGER, INTENT(IN)  :: nk, jspin
      INTEGER, OPTIONAL, INTENT(IN)  :: list(:)
      REAL, INTENT(OUT) :: z(:, :)

      INTEGER :: nmat
      INTEGER i

      nmat = SIZE(z, 1)
      !read eigenvectors
      IF (.NOT. PRESENT(list)) THEN
         ! read all eigenvectors
         CALL io_read_real2(d%evsetid, (/1, 1, 1, nk, jspin/), &
                            (/1, nmat, SIZE(z, 2), 1, 1/), z(:nmat, :))
      ELSE
         DO i = 1, SIZE(list)
            CALL io_read_real1(d%evsetid, (/1, 1, list(i), nk, jspin/),&
                 &                      (/1, nmat, 1, 1, 1/), z(:nmat, i))
         ENDDO
      END IF

   END SUBROUTINE priv_r_vec

#endif

   SUBROUTINE write_eig(id, nk, jspin, neig, neig_total, eig, n_size, n_rank, zmat, smat)

      !*****************************************************************
      !     writes all eignevecs for the nk-th kpoint
      !*****************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN)          :: id, nk, jspin
      INTEGER, INTENT(IN), OPTIONAL :: n_size, n_rank
      INTEGER, INTENT(IN), OPTIONAL :: neig, neig_total
      REAL, INTENT(IN), OPTIONAL :: eig(:)
      TYPE(t_mat), INTENT(IN), OPTIONAL :: zmat, smat

      INTEGER i, j, k, nv_local, n1, n2, ne
      TYPE(t_data_HDF), POINTER::d
      CALL priv_find_data(id, d)

      if(present(smat)) call juDFT_error("writing smat in HDF not supported yet")
#ifdef CPP_HDF
      !
      !write enparas
      !
      nv_local = HUGE(1)

      !
      !write eigenvalues
      !

      IF (PRESENT(neig_total)) THEN
         CALL io_write_integer0(d%neigsetid, (/nk, jspin/), (/1, 1/), neig_total)
      ENDIF

      IF (PRESENT(n_rank) .AND. PRESENT(n_size) .AND.&
           &        PRESENT(eig) .AND. PRESENT(neig)) THEN
         CALL io_write_real1s(&
              &                     d%energysetid, (/n_rank + 1, nk, jspin/),        &
              &                     (/neig, 1, 1/), eig(:neig), (/n_size, 1, 1/))
         !write eigenvectors
         !
      ELSEIF (PRESENT(eig) .AND. PRESENT(neig)) THEN
         CALL io_write_real1s(&
              &                     d%energysetid, (/1, nk, jspin/),&
              &                     (/neig, 1, 1/), eig(:neig), (/1, 1, 1/))
      ELSE
         IF (PRESENT(eig)) CALL juDFT_error("BUG in calling write_eig")
      ENDIF
      IF (PRESENT(zmat) .AND. .NOT. PRESENT(neig))&
           &    CALL juDFT_error("BUG in calling write_eig with eigenvector")

      n1 = 1; n2 = 0
      IF (PRESENT(n_size)) n1 = n_size
      IF (PRESENT(n_rank)) n2 = n_rank
      IF (PRESENT(zmat)) THEN
         IF (zmat%l_real) THEN
            CALL io_write_real2s(&
                 &                     d%evsetid, (/1, 1, n2 + 1, nk, jspin/),&
                 &           (/1, SIZE(zmat%data_r, 1), neig, 1, 1/), REAL(zmat%data_r(:, :neig)), (/1, 1, n1, 1, 1/))
         ELSE
            CALL io_write_real2s(&
                 &                     d%evsetid, (/1, 1, n2 + 1, nk, jspin/),&
                 &           (/1, SIZE(zmat%data_c, 1), neig, 1, 1/), REAL(zmat%data_c(:, :neig)), (/1, 1, n1, 1, 1/))
            CALL io_write_real2s(&
                 &                     d%evsetid, (/2, 1, n2 + 1, nk, jspin/),&
                 &           (/1, SIZE(zmat%data_c, 1), neig, 1, 1/), AIMAG(zmat%data_c(:, :neig)),&
                 &           (/1, 1, n1, 1, 1/))
         ENDIF
      ENDIF

#endif
   END SUBROUTINE write_eig

#ifdef CPP_HDF

   !----------------------------------------------------------------------
   SUBROUTINE priv_r_vecc(&
        &                     d, nk, jspin, list, z)

      USE m_hdf_tools
      IMPLICIT NONE
      TYPE(t_data_HDF), INTENT(IN)::d
      INTEGER, INTENT(IN)  :: nk, jspin
      INTEGER, OPTIONAL, INTENT(IN)  :: list(:)
      COMPLEX, INTENT(OUT) :: z(:, :)

      REAL, ALLOCATABLE :: z1(:, :, :)
      INTEGER i, j
      INTEGER :: nmat

      nmat = SIZE(z, 1)

      IF (.NOT. PRESENT(list)) THEN
         ! read all eigenvectors
         ALLOCATE (z1(2, nmat, SIZE(z, 2)))
         CALL io_read_real3(d%evsetid, (/1, 1, 1, nk, jspin/),&
              &                      (/2, nmat, SIZE(z, 2), 1, 1/), z1)
         DO i = 1, SIZE(z, 2)
            z(:, i) = CMPLX(z1(1, :, i), z1(2, :, i))
         ENDDO
      ELSE
         ALLOCATE (z1(2, nmat, 1))
         DO i = 1, SIZE(list)
            CALL io_read_real3(d%evsetid, (/1, 1, list(i), nk, jspin/),&
             &                      (/2, nmat, 1, 1, 1/), z1)
            z(:, i) = CMPLX(z1(1, :, 1), z1(2, :, 1))
         ENDDO
      END IF
   END SUBROUTINE priv_r_vecc
   !-----------------------------------------------------------------------

#endif

   SUBROUTINE read_eig(id, nk, jspin, neig, eig, list, zMat, smat)
      IMPLICIT NONE
      INTEGER, INTENT(IN)            :: id, nk, jspin
      INTEGER, INTENT(OUT), OPTIONAL  :: neig
      REAL, INTENT(OUT), OPTIONAL  :: eig(:)
      INTEGER, INTENT(IN), OPTIONAL   :: list(:)
      TYPE(t_mat), OPTIONAL  :: zmat, smat

#ifdef CPP_HDF
      INTEGER:: n1, n, k
      TYPE(t_data_HDF), POINTER::d
      if(present(smat)) call juDFT_error("reading smat not supported for HDF")
      CALL priv_find_data(id, d)

      IF (PRESENT(neig)) THEN
         CALL io_read_integer0(d%neigsetid, (/nk, jspin/), (/1, 1/), neig)

         IF (PRESENT(eig)) THEN                           ! read eigenv
            IF (neig > SIZE(eig)) THEN
               WRITE (*, *) neig, SIZE(eig)
               CALL juDFT_error("eig66_hdf$readeig", calledby="eig66_hdf")
            ENDIF
            CALL io_read_real1(d%energysetid, (/1, nk, jspin/), (/neig, 1, 1/),&
                 &                      eig(:neig))
         ENDIF
      ENDIF

      IF (PRESENT(zMat)) THEN
         IF (zmat%l_real) THEN
            CALL priv_r_vec(d, nk, jspin, list, zmat%data_r)
         ELSE
            CALL priv_r_vecc(d, nk, jspin, list, zmat%data_c)
         ENDIF
      ENDIF
#endif
   END SUBROUTINE read_eig

END MODULE

