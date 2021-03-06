!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_mpimat
   USE m_judft
   USE m_types_mat
   USE m_constants
   IMPLICIT NONE
   PRIVATE
   INTEGER, PARAMETER    :: DEFAULT_BLOCKSIZE = 64
   INTEGER, PARAMETER   :: dlen_ = 9

   !<This data-type extends the basic t_mat for distributed matrices.
   !<
   !<It stores the additional mpi_communicator and sets up a blacs grid for the matrix.
   !<This can be used to perform scalapack calls on the matrix with little additional input.
   !<The copy procedure is overwritten from t_mat to enable also redistribution of the matrix.
   TYPE t_blacsdata
      INTEGER:: no_use
      INTEGER:: mpi_com                          !> mpi-communiator over which matrix is distributed
      INTEGER:: blacs_desc(dlen_)                !> blacs descriptor
      !> 1: =1
      !> 2: context
      !> 3,4: global matrix size
      !> 5,6: block sizes
      !> 7,8: row/colum of grid for first row/colum of matrix
      !> 9: leading dimension of local matrix
      INTEGER:: npcol, nprow                     !> the number of columns/rows in the processor grid
      INTEGER:: mycol, myrow
   END TYPE t_blacsdata

   TYPE, EXTENDS(t_mat):: t_mpimat
      INTEGER                   :: global_size1, global_size2        !> this is the size of the full-matrix
      TYPE(t_blacsdata), POINTER :: blacsdata
   CONTAINS
      PROCEDURE   :: copy => mpimat_copy     !<overwriten from t_mat, also performs redistribution
      PROCEDURE   :: move => mpimat_move     !<overwriten from t_mat, also performs redistribution
      PROCEDURE   :: free => mpimat_free     !<overwriten from t_mat, takes care of blacs-grids
      PROCEDURE   :: multiply => mpimat_multiply  !<overwriten from t_mat, takes care of blacs-grids
      PROCEDURE   :: init_details => mpimat_init
      PROCEDURE   :: init_template => mpimat_init_template     !<overwriten from t_mat, also calls alloc in t_mat
      PROCEDURE   :: add_transpose => mpimat_add_transpose !<overwriten from t_mat
      PROCEDURE   :: u2l => t_mpimat_u2l   ! construct full matrix if only upper triangle of hermitian matrix is given
      PROCEDURE   :: l2u => t_mpimat_l2u
      PROCEDURE   :: print_matrix
      PROCEDURE   :: from_non_dist
      PROCEDURE   :: transpose => mpimat_transpose
      procedure   :: print_type => mpimat_print_type
      FINAL :: finalize, finalize_1d, finalize_2d, finalize_3d
   END TYPE t_mpimat

   PUBLIC t_mpimat

CONTAINS
   subroutine mpimat_print_type(mat)
      implicit none 
      CLASS(t_mpimat), INTENT(IN)     :: mat 

      write (*,*) "type -> t_ mpimat"
   end subroutine mpimat_print_type

   SUBROUTINE mpimat_multiply(mat1, mat2, res, transA, transB)
      use m_judft
      CLASS(t_mpimat), INTENT(INOUT)     :: mat1
      CLASS(t_mat), INTENT(IN)           :: mat2
      CLASS(t_mat), INTENT(INOUT), OPTIONAL :: res
      character(len=1), intent(in), optional :: transA, transB

#ifdef CPP_SCALAPACK
      TYPE(t_mpimat)::m, r
      character(len=1)  :: transA_i, transB_i

      transA_i = "N"
      if (present(transA)) transA_i = transA
      transB_i = "N"
      if (present(transB)) transB_i = transB
      if (transA /= "N" .or. transB /= "N") call judft_error("trans /= 'N' not yet implemented for MPI")

      IF (.NOT. PRESENT(res)) CALL judft_error("BUG: in mpicase the multiply requires the optional result argument")
      SELECT TYPE (mat2)
      TYPE IS (t_mpimat)
         SELECT TYPE (res)
         TYPE is (t_mpimat)
            CALL m%init(mat1, mat2%global_size1, mat2%global_size2)
            CALL m%copy(mat2, 1, 1)
            CALL r%init(mat1, res%global_size1, res%global_size2)
            IF (mat1%l_real) THEN
               CALL pdgemm(transA_i, transB_i, mat1%global_size1, m%global_size2, mat1%global_size2, 1.0, &
                           mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, &
                           m%data_r, 1, 1, m%blacsdata%blacs_desc, 0.0, &
                           r%data_r, 1, 1, r%blacsdata%blacs_desc)
            ELSE
               CALL pzgemm(transA_i, transB_i, mat1%global_size1, m%global_size2, mat1%global_size2, cmplx_1, &
                           mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, &
                           m%data_c, 1, 1, m%blacsdata%blacs_desc, cmplx_0, &
                           r%data_c, 1, 1, r%blacsdata%blacs_desc)
            ENDIF
            CALL res%copy(r, 1, 1)
            CALL r%free()
            CALL m%free()
         CLASS default
            CALL judft_error("BUG in mpimat%multiply: res needs to be t_mpimat")
         END SELECT
      CLASS default
         CALL judft_error("BUG in mpimat%multiply: mat2 needs to be t_mpimat")
      END SELECT
#endif
   END SUBROUTINE mpimat_multiply

   subroutine mpimat_transpose(mat1, res)
      CLASS(t_mpimat), INTENT(INOUT) ::mat1
      CLASS(t_mat), INTENT(OUT), OPTIONAL ::res
      real, allocatable :: rd(:, :)
      complex, allocatable :: cd(:, :)
#ifdef CPP_SCALAPACK
      if (present(res)) Then
         select type (res)
         type is (t_mpimat)
            res%blacsdata = mat1%blacsdata
            res%matsize1 = mat1%matsize1
            res%matsize2 = mat1%matsize2
            res%global_size1 = mat1%global_size1
            res%global_size2 = mat1%global_size2
            res%l_real = mat1%l_real
         class default
            call judft_error("BUG in t_mpimat%transpose, wrong matrix type")
         end select
      ENDIF

      IF (mat1%l_real) THEN
         allocate (rd(size(mat1%data_r, 1), size(mat1%data_r, 2)))
         call pdtran(mat1%global_size1, mat1%global_size2, 1.0, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, 0.0, rd, 1, 1, mat1%blacsdata%blacs_desc)
         if (present(res)) Then
            call move_alloc(rd, res%data_r)
         else
            call move_alloc(rd, mat1%data_r)
         endif
      ELSE
         allocate (cd(size(mat1%data_c, 1), size(mat1%data_c, 2)))
         call pztranc(mat1%global_size1, mat1%global_size2, cmplx(1.0, 0.0), mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, cmplx(0.0, 0.0), cd, 1, 1, mat1%blacsdata%blacs_desc)
         if (present(res)) Then
            call move_alloc(cd, res%data_c)
         else
            call move_alloc(cd, mat1%data_c)
         endif
      ENDIF
#else
      call judft_error("Not compiled with MPI")
#endif
   END subroutine

   SUBROUTINE print_matrix(mat, fileno)
      CLASS(t_mpimat), INTENT(INOUT) ::mat
      INTEGER:: fileno

#ifdef CPP_SCALAPACK
      INCLUDE 'mpif.h'
      INTEGER, EXTERNAL:: indxl2g
      CHARACTER(len=10)::filename
      INTEGER :: irank, isize, i, j, npr, npc, r, c, tmp, err, status(MPI_STATUS_SIZE)

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, irank, err)
      CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, isize, err)

      tmp = 0

      IF (irank > 0) CALL MPI_RECV(tmp, 1, MPI_INTEGER, irank - 1, 0, mat%blacsdata%mpi_com, status, err) !lock
      WRITE (filename, "(a,i0)") "out.", fileno
      OPEN (fileno, file=filename, access='append')

      CALL blacs_gridinfo(mat%blacsdata%blacs_desc(2), npr, npc, r, c)
      DO i = 1, mat%matsize1
         DO j = 1, mat%matsize2
            IF (mat%l_real) THEN
               WRITE (fileno, "(5(i0,1x),2(f10.5,1x))") irank, i, j, indxl2g(i, mat%blacsdata%blacs_desc(5), r, 0, npr), &
                  indxl2g(j, mat%blacsdata%blacs_desc(6), c, 0, npc), mat%data_r(i, j)
            ELSE
               WRITE (fileno, "(5(i0,1x),2(f10.5,1x))") irank, i, j, indxl2g(i, mat%blacsdata%blacs_desc(5), r, 0, npr), &
                  indxl2g(j, mat%blacsdata%blacs_desc(6), c, 0, npc), mat%data_c(i, j)
            END IF
         ENDDO
      ENDDO
      CLOSE (fileno)
      IF (irank + 1 < isize) CALL MPI_SEND(tmp, 1, MPI_INTEGER, irank + 1, 0, mat%blacsdata%mpi_com, err)

#endif
   END SUBROUTINE print_matrix

   subroutine t_mpimat_l2u(mat)
      implicit none
      class(t_mpimat), intent(inout) :: mat

      call judft_error("l2u not yet implemented for t_mpimat")
   end subroutine t_mpimat_l2u

   SUBROUTINE t_mpimat_u2l(mat)
      implicit none
      CLASS(t_mpimat), INTENT(INOUT) ::mat

      INTEGER :: i, j, i_glob, j_glob, myid, err, np
      COMPLEX, ALLOCATABLE:: tmp_c(:, :)
      REAL, ALLOCATABLE   :: tmp_r(:, :)
#ifdef CPP_SCALAPACK
      INCLUDE 'mpif.h'
      INTEGER, EXTERNAL    :: numroc, indxl2g  !SCALAPACK functions

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, myid, err)
      CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, np, err)
      !Set lower part of matrix to zero

      DO i = 1, mat%matsize1
         DO j = 1, mat%matsize2
            ! Get global column corresponding to i and number of local rows up to
            ! and including the diagonal, these are unchanged in A
            i_glob = indxl2g(i, mat%blacsdata%blacs_desc(5), mat%blacsdata%myrow, 0, mat%blacsdata%nprow)
            j_glob = indxl2g(j, mat%blacsdata%blacs_desc(6), mat%blacsdata%mycol, 0, mat%blacsdata%npcol)

            IF (i_glob > j_glob) THEN
               IF (mat%l_real) THEN
                  mat%data_r(i, j) = 0.0
               ELSE
                  mat%data_c(i, j) = 0.0
               ENDIF
            ENDIF
            IF (i_glob == j_glob) THEN
               IF (mat%l_real) THEN
                  mat%data_r(i, j) = mat%data_r(i, j)/2.0
               ELSE
                  mat%data_c(i, j) = mat%data_c(i, j)/2.0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (mat%l_real) THEN
         ALLOCATE (tmp_r(mat%matsize1, mat%matsize2))
         tmp_r = mat%data_r
      ELSE
         ALLOCATE (tmp_c(mat%matsize1, mat%matsize2))
         tmp_c = mat%data_c
      END IF
      CALL MPI_BARRIER(mat%blacsdata%mpi_com, i)
      IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK

         CALL pdgeadd('t', mat%global_size1, mat%global_size2, 1.0, tmp_r, 1, 1, mat%blacsdata%blacs_desc, 1.0, mat%data_r, 1, 1, mat%blacsdata%blacs_desc)
      ELSE
         CALL pzgeadd('c', mat%global_size1, mat%global_size2, CMPLX(1.0, 0.0), tmp_c, 1, 1, mat%blacsdata%blacs_desc, CMPLX(1.0, 0.0), mat%data_c, 1, 1, mat%blacsdata%blacs_desc)
#endif
      END IF
#endif
   END SUBROUTINE t_mpimat_u2l

   SUBROUTINE mpimat_add_transpose(mat, mat1)
      CLASS(t_mpimat), INTENT(INOUT) ::mat
      CLASS(t_mat), INTENT(INOUT) ::mat1

      INTEGER:: i, ii, n_size, n_rank

      SELECT TYPE (mat1)
      TYPE IS (t_mpimat)
#ifdef CPP_MPI
         CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, n_rank, i)
         CALL MPI_COMM_SIZE(mat%blacsdata%mpi_com, n_size, i)
#endif
         !Set lower part of matrix to zero...
         ii = 0
         DO i = n_rank + 1, MIN(mat%global_size1, mat%global_size2), n_size
            ii = ii + 1
            IF (mat%l_real) THEN
               mat%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i + 1:, ii) = 0.0
               mat1%data_r(i, ii) = 0.0
            ELSE
               mat%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i + 1:, ii) = 0.0
               mat1%data_c(i, ii) = 0.0
            ENDIF
         ENDDO
         IF (mat%l_real) THEN
#ifdef CPP_SCALAPACK

            CALL pdgeadd('t', mat1%global_size1, mat1%global_size2, 1.0, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, 1.0, mat%data_r, 1, 1, mat%blacsdata%blacs_desc)
         ELSE
            CALL pzgeadd('c', mat1%global_size1, mat1%global_size2, CMPLX(1.0, 0.0), mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, CMPLX(1.0, 0.0), mat%data_c, 1, 1, mat1%blacsdata%blacs_desc)
#endif
         END IF
         !Now multiply the diagonal of the matrix by 1/2

         !ii=0
         !DO i=n_rank+1,MIN(mat%global_size1,mat%global_size2),n_size
         !   ii=ii+1
         !   IF (mat%l_real) THEN
         !      mat%data_r(i,ii)=mat%data_r(i,ii)/2
         !   ELSE
         !      mat%data_c(i,ii)=mat%data_c(i,ii)/2
         !   END IF
         !ENDDO
      CLASS default
         CALL judft_error("Inconsistent types in t_mpimat_add_transpose")
      END SELECT

   END SUBROUTINE mpimat_add_transpose

   SUBROUTINE mpimat_copy(mat, mat1, n1, n2)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      CLASS(t_mat), INTENT(IN)      ::mat1
      INTEGER, INTENT(IN) ::n1, n2

      select type (mat1)
      type is(t_mpimat)
          
      class default
         call judft_error("you can only copy a t_mpimat to a t_mpimat")
      end select

#ifdef CPP_SCALAPACK
      SELECT TYPE (mat1)
      TYPE IS (t_mpimat)
         IF (mat%l_real) THEN
            CALL pdgemr2d(Mat1%global_size1, mat1%global_size2, mat1%data_r, 1, 1, mat1%blacsdata%blacs_desc, mat%data_r, n1, n2, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
         ELSE
            CALL pzgemr2d(mat1%global_size1, mat1%global_size2, mat1%data_c, 1, 1, mat1%blacsdata%blacs_desc, mat%data_c, n1, n2, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
         END IF
      CLASS DEFAULT
         CALL judft_error("Wrong datatype in copy")
      END SELECT
#endif
   END SUBROUTINE mpimat_copy

   SUBROUTINE from_non_dist(mat, mat1)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      TYPE(t_mat), INTENT(IN)       ::mat1

      INTEGER:: blacs_desc(9), irank, ierr, umap(1, 1), np
#ifdef CPP_SCALAPACK
      blacs_desc = (/1, -1, mat1%matsize1, mat1%matsize2, mat1%matsize1, mat1%matsize2, 0, 0, mat1%matsize1/)

      CALL MPI_COMM_RANK(mat%blacsdata%mpi_com, irank, ierr)
      umap(1, 1) = 0
      CALL BLACS_GET(mat%blacsdata%blacs_desc(2), 10, blacs_desc(2))
      CALL BLACS_GRIDMAP(blacs_desc(2), umap, 1, 1, 1)
      IF (mat%l_real) THEN
         CALL pdgemr2d(Mat1%matsize1, mat1%matsize2, mat1%data_r, 1, 1, blacs_desc, mat%data_r, 1, 1, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
      ELSE
         CALL pzgemr2d(mat1%matsize1, mat1%matsize2, mat1%data_c, 1, 1, blacs_desc, mat%data_c, 1, 1, mat%blacsdata%blacs_desc, mat%blacsdata%blacs_desc(2))
      END IF
#endif
   END SUBROUTINE from_non_dist

   SUBROUTINE mpimat_move(mat, mat1)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)::mat
      CLASS(t_mat), INTENT(INOUT)   ::mat1
      CALL mat%copy(mat1, 1, 1)
   END SUBROUTINE mpimat_move

   SUBROUTINE finalize(mat)
      IMPLICIT NONE
      TYPE(t_mpimat), INTENT(INOUT) :: mat
      CALL mpimat_free(mat)
   END SUBROUTINE finalize

   SUBROUTINE finalize_1d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:)
      INTEGER                      :: i
      DO i = 1, size(mat)
         CALL mpimat_free(mat(i))
      ENDDO
   END SUBROUTINE finalize_1d

   SUBROUTINE finalize_2d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:, :)
      INTEGER                      :: i, j

      DO i = 1, size(mat, dim=1)
         DO j = 1, size(mat, dim=2)
            CALL mpimat_free(mat(i, j))
         ENDDO
      ENDDO
   END SUBROUTINE finalize_2d

   SUBROUTINE finalize_3d(mat)
      IMPLICIT NONE

      TYPE(t_mpimat), INTENT(INOUT) :: mat(:, :, :)
      INTEGER                      :: i, j, k

      DO i = 1, size(mat, dim=1)
         DO j = 1, size(mat, dim=2)
            DO k = 1, size(mat, dim=3)
               CALL mpimat_free(mat(i, j, k))
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE finalize_3d

   SUBROUTINE mpimat_free(mat)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT) :: mat
      INTEGER :: ierr
      IF (ALLOCATED(mat%data_r)) DEALLOCATE (mat%data_r)
      IF (ALLOCATED(mat%data_c)) DEALLOCATE (mat%data_c)
      IF (ASSOCIATED(mat%blacsdata)) THEN
         IF (mat%blacsdata%no_use > 1) THEN
            mat%blacsdata%no_use = mat%blacsdata%no_use - 1
            mat%blacsdata => null()
         ELSE
#ifdef CPP_SCALAPACK
            CALL BLACS_GRIDEXIT(mat%blacsdata%blacs_desc(2), ierr)
            DEALLOCATE (mat%blacsdata)
#endif
         END IF
      ENDIF
   END SUBROUTINE mpimat_free

   !>Initialization of the distributed matrix.
  !!
  !! The argument l_2d controls the kind of distribution used:
  !!  - TRUE: the matrix is a Scalapack BLOCK-CYCLIC distribution
  !!  - FALSE: the matrix is distributed in a one-dimensional column cyclic distribution with blocksize 1
  !! as used in the parallel matrix setup of FLEUR
   SUBROUTINE mpimat_init(mat, l_real, matsize1, matsize2, mpi_subcom, l_2d, nb_x, nb_y)
      IMPLICIT NONE
      CLASS(t_mpimat)             :: mat
      INTEGER, INTENT(IN), OPTIONAL :: matsize1, matsize2, mpi_subcom
      LOGICAL, INTENT(IN), OPTIONAL :: l_real, l_2d
      INTEGER, INTENT(IN), OPTIONAL :: nb_y, nb_x
#ifdef CPP_SCALAPACK
      INTEGER::nbx, nby, irank, ierr
      include 'mpif.h'
      nbx = DEFAULT_BLOCKSIZE; nby = DEFAULT_BLOCKSIZE
      IF (PRESENT(nb_x)) nbx = nb_x
      IF (PRESENT(nb_y)) nby = nb_y
      IF (.NOT. (PRESENT(matsize1) .AND. PRESENT(matsize2) .AND. PRESENT(mpi_subcom) .AND. PRESENT(l_real) .AND. PRESENT(l_2d))) &
         CALL judft_error("Optional arguments must be present in mpimat_init")
      mat%global_size1 = matsize1
      mat%global_size2 = matsize2
      ALLOCATE (mat%blacsdata)
      mat%blacsdata%no_use = 1
      CALL priv_create_blacsgrid(mpi_subcom, l_2d, matsize1, matsize2, nbx, nby, &
                                 mat%blacsdata, mat%matsize1, mat%matsize2)
      mat%blacsdata%mpi_com = mpi_subcom
      CALL mat%alloc(l_real) !Attention,sizes determined in call to priv_create_blacsgrid
      !check if this matrix is actually distributed over MPI_COMM_SELF
      IF (mpi_subcom == MPI_COMM_SELF) THEN
         CALL MPI_COMM_RANK(mpi_subcom, irank, ierr)
         IF (irank > 0) mat%blacsdata%blacs_desc(2) = -1
      END IF
#endif
   END SUBROUTINE mpimat_init

   SUBROUTINE mpimat_init_template(mat, templ, global_size1, global_size2)
      IMPLICIT NONE
      CLASS(t_mpimat), INTENT(INOUT)  :: mat
      CLASS(t_mat), INTENT(IN)        :: templ
      INTEGER, INTENT(IN), OPTIONAL    :: global_size1, global_size2

      INTEGER::numroc
      EXTERNAL::numroc

      SELECT TYPE (templ)
      TYPE IS (t_mpimat)
         mat%l_real = templ%l_real
         IF (PRESENT(global_size1) .AND. PRESENT(global_size2)) THEN
            ALLOCATE (mat%blacsdata)
            mat%blacsdata = templ%blacsdata
            mat%blacsdata%no_use = 1
            mat%blacsdata%blacs_desc(3) = global_size1
            mat%blacsdata%blacs_desc(4) = global_size2
            mat%global_size1 = global_size1
            mat%global_size2 = global_size2
#ifdef CPP_SCALAPACK
            mat%matsize1 = NUMROC(global_size1, mat%blacsdata%blacs_desc(5), mat%blacsdata%myrow, mat%blacsdata%blacs_desc(7), mat%blacsdata%nprow)
            mat%matsize1 = NUMROC(global_size2, mat%blacsdata%blacs_desc(6), mat%blacsdata%mycol, mat%blacsdata%blacs_desc(8), mat%blacsdata%npcol)
#endif
         ELSE
            mat%matsize1 = templ%matsize1
            mat%matsize2 = templ%matsize2
            mat%global_size1 = templ%global_size1
            mat%global_size2 = templ%global_size2
            mat%blacsdata => templ%blacsdata
            mat%blacsdata%no_use = mat%blacsdata%no_use + 1
         ENDIF
         CALL mat%alloc()

      CLASS default
         CALL judft_error("Mixed initialization in t_mpimat not possible(BUG)")
      END SELECT
   END SUBROUTINE mpimat_init_template

   SUBROUTINE priv_create_blacsgrid(mpi_subcom, l_2d, m1, m2, nbc, nbr, blacsdata, local_size1, local_size2)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mpi_subcom
      INTEGER, INTENT(IN) :: m1, m2
      INTEGER, INTENT(INOUT)::nbc, nbr
      LOGICAL, INTENT(IN) :: l_2d
      type(t_blacsdata), INTENT(OUT)::blacsdata
      INTEGER, INTENT(OUT):: local_size1, local_size2

#ifdef CPP_SCALAPACK
      INCLUDE 'mpif.h'
      INTEGER     :: myrowssca, mycolssca
      INTEGER     :: iamblacs, npblacs, np, myid, mycol, myrow
      INTEGER     :: nprow2, npcol2
      INTEGER     :: k, i, j, ictextblacs
      INTEGER     :: ierr

      INTEGER, ALLOCATABLE :: iblacsnums(:), ihelp(:), iusermap(:, :)

      EXTERNAL descinit, blacs_get
      EXTERNAL blacs_pinfo, blacs_gridinit

      !Determine rank and no of processors
      CALL MPI_COMM_RANK(mpi_subcom, myid, ierr)
      CALL MPI_COMM_SIZE(mpi_subcom, np, ierr)

      ! compute processor grid, as square as possible
      ! If not square with more rows than columns
      IF (l_2d) THEN
         distloop: DO j = INT(SQRT(REAL(np))), 1, -1
            IF ((np/j)*j == np) THEN
               blacsdata%npcol = np/j
               blacsdata%nprow = j
               EXIT distloop
            ENDIF
         ENDDO distloop
      ELSE
         nbc = 1
         nbr = MAX(m1, m2)
         blacsdata%npcol = np
         blacsdata%nprow = 1
      ENDIF
      ALLOCATE (iblacsnums(np), ihelp(np), iusermap(blacsdata%nprow, blacsdata%npcol))

      !   A blacsdata%nprow*blacsdata%npcol processor grid will be created
      !   Row and column index myrow, mycol of this processor in the grid
      !   and distribution of A and B in ScaLAPACK
      !   The local processor will get myrowssca rows and mycolssca columns
      !   of A and B
      !

      myrow = myid/blacsdata%npcol  ! my row number in the BLACS blacsdata%nprow*blacsdata%npcol grid
      mycol = myid - (myid/blacsdata%npcol)*blacsdata%npcol  ! my column number in the BLACS blacsdata%nprow*blacsdata%npcol grid
      !
      !  Now allocate Asca to put the elements of Achi or receivebuffer to
      !
      myrowssca = (m1 - 1)/(nbr*blacsdata%nprow)*nbr + MIN(MAX(m1 - (m1 - 1)/(nbr*blacsdata%nprow)*nbr*blacsdata%nprow - nbr*myrow, 0), nbr)
      !     Number of rows the local process gets in ScaLAPACK distribution
      mycolssca = (m2 - 1)/(nbc*blacsdata%npcol)*nbc + MIN(MAX(m2 - (m2 - 1)/(nbc*blacsdata%npcol)*nbc*blacsdata%npcol - nbc*mycol, 0), nbc)

      !Get BLACS ranks for all MPI ranks
      CALL BLACS_PINFO(iamblacs, npblacs)  ! iamblacs = local process rank (e.g. myid)
      ! npblacs  = number of available processes
      iblacsnums = -2
      ihelp = -2
      ihelp(myid + 1) = iamblacs ! Get the Blacs id corresponding to the MPI id
      CALL MPI_ALLREDUCE(ihelp, iblacsnums, np, MPI_INTEGER, MPI_MAX, mpi_subcom, ierr)
      IF (ierr /= 0) call judft_error('Error in allreduce for BLACS nums')

      !     iblacsnums(i) is the BLACS-process number of MPI-process i-1
      k = 1
      DO i = 1, blacsdata%nprow
         DO j = 1, blacsdata%npcol
            iusermap(i, j) = iblacsnums(k)
            k = k + 1
         ENDDO
      ENDDO
!#ifdef CPP_BLACSDEFAULT
      !Get the Blacs default context
      CALL BLACS_GET(0, 0, ictextblacs)
!#else
!    ictextblacs=mpi_subcom
!#endif
      ! Create the Grid
      CALL BLACS_GRIDMAP(ictextblacs, iusermap, size(iusermap, 1), blacsdata%nprow, blacsdata%npcol)
      !     Now control, whether the BLACS grid is the one we wanted
      CALL BLACS_GRIDINFO(ictextblacs, nprow2, npcol2, blacsdata%myrow, blacsdata%mycol)
      IF (nprow2 /= blacsdata%nprow) THEN
         WRITE (oUnit, *) 'Wrong number of rows in BLACS grid'
         WRITE (oUnit, *) 'nprow=', blacsdata%nprow, ' nprow2=', nprow2
         call judft_error('Wrong number of rows in BLACS grid')
      ENDIF
      IF (npcol2 /= blacsdata%npcol) THEN
         WRITE (oUnit, *) 'Wrong number of columns in BLACS grid'
         WRITE (oUnit, *) 'npcol=', blacsdata%npcol, ' npcol2=', npcol2
         call judft_error('Wrong number of columns in BLACS grid')

      ENDIF

      !Create the descriptors
      CALL descinit(blacsdata%blacs_desc, m1, m2, nbr, nbc, 0, 0, ictextblacs, myrowssca, ierr)
      IF (ierr /= 0) call judft_error('Creation of BLACS descriptor failed')
      local_size1 = myrowssca
      local_size2 = mycolssca
#endif
   END SUBROUTINE priv_create_blacsgrid
END MODULE m_types_mpimat
