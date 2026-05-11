!
!     isop4_axis95  Elementos finitos elasticos axissimetricos
!     ============  isoparametricos quadrilateros de quatro pontos nodais
!
!     Adaptado de:
!     BREBBIA,C.A. & FERRANTE,A.J. (1986) Computational methods for the
!     solution of engineering problems. London: Pentech Press. 370p.
!
!     SMITH,I.M. & GRIFFITHS,D.V. (1988) Programming the finite element
!     method. Second Edition. Chichester: John Wiley & Sons. 469p.
!
!     P¢s-processamento no GiD https://www.gidhome.com/
!
module variaveis
implicit none

  integer, parameter :: ndf=2 !! numero de graus de liberdade por n¢
  integer, parameter :: nne=4 !! numero de n¢s por elemento
  integer, parameter :: ndfel=ndf*nne !! graus de liberdade por elemento

  integer, parameter :: in=15 !! numero do arquivo de entrada
  integer, parameter :: io=16 !! numero do arquivo de saida
  integer, parameter :: i33=33 !! numero do arquivo do GiD
  integer, parameter :: i34=34 !! numero do arquivo do GiD

  integer :: nnode !! numero de n¢s
  real(4), allocatable :: X(:),Y(:) !! coordenadas dos n¢s
  integer :: nelem !! numero de elementos
  integer, allocatable :: kon(:) !! conectividade dos elementos
  integer, allocatable :: imat(:) !! numero do material do elemento
  integer :: nmat !! numero de materiais
  real(4), allocatable :: Ei(:) !! modulo de elasticidade do elemento
  real(4), allocatable :: Poisson(:) !! Poisson do elemento
  integer :: nln !! numero de nos carregados
  integer :: nbn !! numero de nos com condicao de contorno
  integer, allocatable :: istatus(:) !! indicador de status
  real(4), allocatable :: Prescrito(:) !! deslocamentos prescritos

  integer :: neq !! numero total de incognitas
  real(4), allocatable :: Carga(:) !! Vetor de cargas nodais
  real(4), allocatable :: Desloc(:) !! Vetor de deslocamentos nodais
  real(4), allocatable :: RigGlobal(:,:) !! Matriz de rigidez global
  real(4) :: RigElem(ndfel,ndfel) !! Matriz de rigidez do elemento

  real(4), allocatable :: stress(:,:,:) !! Tensoes nos pontos de Gauss

end module variaveis



module Entrada_Saida

contains

subroutine Abre_Arquivos()
use variaveis

   open(UNIT=in,file='NomeArq.dat',status='old')
   open(UNIT=io,file='NomeArq.out',status='UNKNOWN')
   open(UNIT=i33,file='NomeArq.post.msh',status='UNKNOWN')
   open(UNIT=i34,file='NomeArq.post.res',status='UNKNOWN')

   return
end subroutine Abre_Arquivos



subroutine Entrada_Dados()
use variaveis
implicit none
integer :: i,j,ic(nne),n1,k,k1,L1,L2,n2
real(4) :: w(ndf)

   write(io,'(A,/)')' DADOS DE ENTRADA'

   read(in,*)nnode
   write(io,'(A,i5)')' numeros de nos          :',nnode
   write(io,'(/,A)')' coordenadas nodais'
   write(io,'(7x,A,6x,A,9x,A)')' no ','x','y'
   allocate(X(nnode)); allocate(Y(nnode))
   do j=1,nnode
      read(in,*)i,X(i),Y(i)
      write(io,'(i10,2f10.2)')i,X(i),Y(i)
   enddo

   read(in,*)nelem
   write(io,'(/,A,i5/)')' numero de elementos     :',nelem
   write(io,*)' conectividade dos elementos e propriedades'
   write(io,'(A,4(7x,A),5x,A)')'elemento','no1','no2','no3','no4','material'
   allocate(kon(nne*nelem))
   allocate(imat(nelem))
   do j=1,nelem
      read(in,*)i,ic(1),ic(2),ic(3),ic(4),imat(i)
      write (io,'(6i10)') i,ic(1),ic(2),ic(3),ic(4),imat(i)
      n1=nne*(i-1)
      kon(n1+1)=ic(1)
      kon(n1+2)=ic(2)
      kon(n1+3)=ic(3)
      kon(n1+4)=ic(4)
   enddo

   read(in,*)nmat
   write(io,'(/,A,i5/)')' numero de materiais     :',nmat
   allocate(Ei(nmat))
   allocate(Poisson(nmat))
   do j=1,nmat
      read(in,*)i,Ei(i),Poisson(i)
      write(io,'(i10,2f20.5)')i,Ei(i),Poisson(i)
   enddo

   neq=nnode*ndf  !! calcula o numero total de incognitas
   allocate(Carga(neq))
   allocate(Desloc(neq))
   do i=1,neq
      Carga(i)=0.0 !! Zera vetor de cargas nodais
   enddo

   read(in,*)nln
   write(io,'(/,A,i5/)')' numero de nos carregados:',nln
   write(io,'(A,/7x,A,7x,A,8x,A)')'cargas nodais','no','px','py'
   do i=1,nln
      read(in,*) j,(w(k),k=1,ndf)
      write(io,'(i10,2f10.2)') j,(w(k),k=1,ndf)
      do k=1,ndf
         k1=ndf*(j-1)+k
         Carga(k1)=w(k)
      enddo
   enddo

   read(in,*)nbn
   write(io,'(//,A,i5//)')' numero de nos suportados:',nbn
   write(io,*)' dados das condicoes de contorno'
   write(io,'(23X,A,14X,A)')'status','valores prescritos'
   write(io,'(16X,A)')'(1:prescrito, 0:livre)'
   write(io,'(7x,A,10x,A,9x,A,16x,A,9x,A)')'no','u','v','u','v'
   allocate(istatus((ndf+1)*nbn))
   allocate(Prescrito(neq))
   do i=1,nbn
      read(in,*) j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      write(io,'(3i10,10x,2f10.4)') j,(ic(k),k=1,ndf),(w(k),k=1,ndf)
      L1=(ndf+1)*(i-1)+1
      L2=ndf*(j-1)
      istatus(L1)=j
      do k=1,ndf
         n1=L1+k
         n2=L2+k
         istatus(n1)=ic(k)
         Prescrito(n2)=w(k)
         !write(io,*)ndf*(nnode-1)+ndf,'L1',L1,'L2',L2,'n1',n1,'n2',n2
      enddo
   enddo

return
end subroutine Entrada_Dados



subroutine Imprime_Resultados()
use variaveis
implicit none
integer :: i,j,k1,k2

      write(io,'(//,A,/)')' RESULTADOS'

      write(io,'(A,/)')' deslocamentos nodais'
      write(io,'(8X,A,12X,A,14x,A)')'no','u','v'
      do i=1,nnode
        k1=ndf*(i-1)+1
        k2=k1+ndf-1
        write(io,'(i10,2E15.5)') i,(Desloc(j),j=k1,k2)
      enddo

      WRITE(io,'(//,A,4X,4(A,8X),/)')'ELEMENT NUMBER','X-STRESS','Y-STRESS','XY-STRESS','ZZ-STRESS'
      do i=1,nelem
        do j=1,4
          write(io,'(6X,I4,5X,4e16.7)')i,stress(i,j,1),stress(i,j,2),stress(i,j,3),stress(i,j,4)
        enddo
      enddo

   return
end subroutine Imprime_Resultados



subroutine Resultados_GiD()
use variaveis
implicit none
integer :: i,j,n1,k1,k2

   !open(UNIT=i33,file=arq1,status='UNKNOWN')
   !! Escreve arquivo da malha
   write(i33,*)'MESH "Malha" dimension 2 ElemType Quadrilateral Nnode 4'
   write(i33,*)'Coordinates'
   do i=1,nnode
      write(i33,*)i,'  ',x(i),'  ',y(i)
   enddo
   write(i33,*)'end coordinates'
   write(i33,*)'Elements'
   do i=1,nelem
      n1=nne*(i-1)
      write(i33,'(6(2X,I8))') i,kon(n1+1),kon(n1+2),kon(n1+3),kon(n1+4),imat(i)
   enddo
   write(i33,*)'end elements'
   !close(UNIT=i33)

   !open(UNIT=i34,file=arq2,status='UNKNOWN')
   !! Escreve arquivo de resultados
   write(i34,*)'GiD Post Results File 1.0'
   write(i34,*)'GaussPoints "Malha_gauss" ElemType Quadrilateral "Malha"'
      write(i34,*)'Number Of Gauss Points: 4'
      write(i34,*)'Natural Coordinates: Internal'
!      write(34,*)'Natural Coordinates: Given'
!      write(34,*)-1./SQRT(3.),'     ',-1./SQRT(3.)
!      write(34,*)+1./SQRT(3.),'     ',-1./SQRT(3.)
!      write(34,*)+1./SQRT(3.),'     ',+1./SQRT(3.)
!      write(34,*)-1./SQRT(3.),'     ',+1./SQRT(3.)
      write(i34,*)'End gausspoints'
   write(i34,*)'Result "Desloc" "Load Analysis" ',1,' Vector OnNodes'
   write(i34,*)'ComponentNames "X-Desloc", "Y-Desloc"'
   write(i34,*)'Values'
   do i=1,nnode
      k1=ndf*(i-1)+1
      k2=k1+ndf-1
      write(i34,'(i10,2(2X,E15.5))')i,(Desloc(j),j=k1,k2)
   enddo
   write(i34,*)'End Values'
      write(i34,*)'Result "Gauss Stress" "Load Analysis" ',1,' ', &
     &'PlainDeformationMatrix OnGaussPoints "Malha_gauss"'
      write(i34,*)'Values'
      do i=1,nelem
        do j=1,4
          if (j.eq.1) then
            write(i34,'(6X,I4,5X,4e16.7)')i,stress(i,j,1),stress(i,j,2),stress(i,j,3),stress(i,j,4)
          else
            write(i34,'(10X,5X,4e16.7)')stress(i,j,1),stress(i,j,2),stress(i,j,4)
          endif
        enddo
      enddo
      write(34,*)'End Values'
   !close(UNIT=i34)

   return
end subroutine Resultados_GiD



subroutine Fecha_Arquivos()
use variaveis

   close(UNIT=in)
   close(UNIT=io)
   close(UNIT=i33)
   close(UNIT=i34)

   return
end subroutine Fecha_Arquivos


end module Entrada_Saida



module Matriz_Elemento

contains

      SUBROUTINE FMDRAD(DEE,IDEE,E,V)
!      THIS SUBROUTINE FORMS THE ELASTIC AXISYMMETRIC
!      STRESS/STRAIN MATRIX
      REAL DEE(IDEE,*)
      V1=1.-V
      C=E/((1.+V)*(1.-2.*V))
      DEE(1,1)=V1*C
      DEE(2,2)=V1*C
      DEE(3,3)=.5*C*(1.-2.*V)
      DEE(4,4)=V1*C
      DEE(1,2)=V*C
      DEE(2,1)=V*C
      DEE(1,3)=0.
      DEE(3,1)=0.
      DEE(1,4)=V*C
      DEE(4,1)=V*C
      DEE(2,3)=0.
      DEE(3,2)=0.
      DEE(2,4)=V*C
      DEE(4,2)=V*C
      DEE(4,3)=0.
      DEE(3,4)=0.
      RETURN
      END SUBROUTINE FMDRAD



      SUBROUTINE GAUSS(SAMP,ISAMP,NGP)
!      THIS SUBROUTINE PROVIDES THE WEIGHTS AND SAMPLING POINTS
!      FOR GAUSS-LEGENDRE QUADRATURE
      REAL SAMP(ISAMP,2)
      select case (NGP)
      case (1)
        SAMP(1,1)=0.
        SAMP(1,2)=2.
      case (2)
        SAMP(1,1)=1./SQRT(3.)
        SAMP(2,1)=-SAMP(1,1)
        SAMP(1,2)=1.
        SAMP(2,2)=1.
      case (3)
        SAMP(1,1)=.2*SQRT(15.)
        SAMP(2,1)=.0
        SAMP(3,1)=-SAMP(1,1)
        SAMP(1,2)=5./9.
        SAMP(2,2)=8./9.
        SAMP(3,2)=SAMP(1,2)
      case (4)
        SAMP(1,1)=.861136311594053
        SAMP(2,1)=.339981043584856
        SAMP(3,1)=-SAMP(2,1)
        SAMP(4,1)=-SAMP(1,1)
        SAMP(1,2)=.347854845137454
        SAMP(2,2)=.652145154862546
        SAMP(3,2)=SAMP(2,2)
        SAMP(4,2)=SAMP(1,2)
      case (5)
        SAMP(1,1)=.906179845938664
        SAMP(2,1)=.538469310105683
        SAMP(3,1)=.0
        SAMP(4,1)=-SAMP(2,1)
        SAMP(5,1)=-SAMP(1,1)
        SAMP(1,2)=.236926885056189
        SAMP(2,2)=.478628670499366
        SAMP(3,2)=.568888888888889
        SAMP(4,2)=SAMP(2,2)
        SAMP(5,2)=SAMP(1,2)
      case (6)
        SAMP(1,1)=.932469514203152
        SAMP(2,1)=.661209386466265
        SAMP(3,1)=.238619186083197
        SAMP(4,1)=-SAMP(3,1)
        SAMP(5,1)=-SAMP(2,1)
        SAMP(6,1)=-SAMP(1,1)
        SAMP(1,2)=.171324492379170
        SAMP(2,2)=.360761573048139
        SAMP(3,2)=.467913934572691
        SAMP(4,2)=SAMP(3,2)
        SAMP(5,2)=SAMP(2,2)
        SAMP(6,2)=SAMP(1,2)
      case (7)
        SAMP(1,1)=.949107912342759
        SAMP(2,1)=.741531185599394
        SAMP(3,1)=.405845151377397
        SAMP(4,1)=.0
        SAMP(5,1)=-SAMP(3,1)
        SAMP(6,1)=-SAMP(2,1)
        SAMP(7,1)=-SAMP(1,1)
        SAMP(1,2)=.129484966168870
        SAMP(2,2)=.279705391489277
        SAMP(3,2)=.381830050505119
        SAMP(4,2)=.417959183673469
        SAMP(5,2)=SAMP(3,2)
        SAMP(6,2)=SAMP(2,2)
        SAMP(7,2)=SAMP(1,2)
      end select
      RETURN
      END SUBROUTINE GAUSS



      SUBROUTINE NULL(A,IA,M,N)
!      THIS SUBROUTINE NULLS A 2-D ARRAY
      REAL A(IA,*)
      do I=1,M
        do J=1,N
          A(I,J)=0.0
        enddo
      enddo
      RETURN
      END SUBROUTINE NULL



      SUBROUTINE FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)
!      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
!      THEIR DERIVATIVES FOR 4-NODED QUADRILATERAL ELEMENTS
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI=SAMP(J,1)
      ETAM=.25*(1.-ETA)
      ETAP=.25*(1.+ETA)
      XIM=.25*(1.-XI)
      XIP=.25*(1.+XI)
      FUN(4)=4.*XIM*ETAM
      FUN(3)=4.*XIM*ETAP
      FUN(2)=4.*XIP*ETAP
      FUN(1)=4.*XIP*ETAM
      DER(1,4)=-ETAM
      DER(1,3)=-ETAP
      DER(1,2)=ETAP
      DER(1,1)=ETAM
      DER(2,4)=-XIM
      DER(2,3)=XIM
      DER(2,2)=XIP
      DER(2,1)=-XIP
      RETURN
      END SUBROUTINE FORMLN



      SUBROUTINE XMATMUL(A,IA,B,IB,C,IC,L,M,N)
!      THIS SUBROUTINE FORMS THE PRODUCT OF TWO MATRICES
      REAL A(IA,*),B(IB,*),C(IC,*)
      do I=1,L
        do J=1,N
          X=0.0
          do K=1,M
            X=X+A(I,K)*B(K,J)
          enddo
          C(I,J)=X
        enddo
      enddo
      RETURN
      END SUBROUTINE XMATMUL



      SUBROUTINE TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
!      THIS SUBROUTINE FORMS THE INVERSE OF A 2 BY 2 MATRIX
      REAL JAC(IJAC,*),JAC1(IJAC1,*)
      DET=JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)
      JAC1(1,1)=JAC(2,2)
      JAC1(1,2)=-JAC(1,2)
      JAC1(2,1)=-JAC(2,1)
      JAC1(2,2)=JAC(1,1)
      do K=1,2
        do L=1,2
          JAC1(K,L)=JAC1(K,L)/DET
        enddo
      enddo
      RETURN
      END SUBROUTINE TWOBY2



      SUBROUTINE FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
!      THIS SUBROUTINE FORMS THE STRAIN/DISPLACEMENT MATRIX
!      FOR AXISYMMETRIC STRAIN
      REAL BEE(IBEE,*),DERIV(IDERIV,*),FUN(*),COORD(ICOORD,*)
      SUM=0.
      do I=1,NOD
        SUM=SUM+FUN(I)*COORD(I,1)
      enddo
      do M=1,NOD
        K=2*M
        L=K-1
        X=DERIV(1,M)
        BEE(1,L)=X
        BEE(3,K)=X
        Y=DERIV(2,M)
        BEE(2,K)=Y
        BEE(3,L)=Y
        BEE(4,L)=FUN(M)/SUM
      enddo
      RETURN
      END SUBROUTINE FMBRAD



      SUBROUTINE MATRAN(A,IA,B,IB,M,N)
!      THIS SUBROUTINE FORMS THE TRANSPOSE OF A MATRIX
      REAL A(IA,*),B(IB,*)
      do I=1,M
        do J=1,N
          A(J,I)=B(I,J)
        enddo
      enddo
      RETURN
      END SUBROUTINE MATRAN



      SUBROUTINE MSMULT(A,IA,C,M,N)
!      THIS SUBROUTINE MULTIPLIES A MATRIX BY A SCALAR
      REAL A(IA,*)
      do I=1,M
        do J=1,N
          A(I,J)=A(I,J)*C
        enddo
      enddo
      RETURN
      END SUBROUTINE MSMULT



      SUBROUTINE MATADD(A,IA,B,IB,M,N)
!      THIS SUBROUTINE ADDS TWO EQUAL SIZED ARRAYS
      REAL A(IA,*),B(IB,*)
      do I=1,M
        do J=1,N
          A(I,J)=A(I,J)+B(I,J)
        enddo
      enddo
      RETURN
      END SUBROUTINE MATADD



      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
!      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
      REAL M(IM,*),V(*),Y(*)
      do I=1,K
        X=0.
        do J=1,L
          X=X+M(I,J)*V(J)
        enddo
        Y(I)=X
      enddo
      RETURN
      END SUBROUTINE MVMULT



subroutine Monta_Matriz_Elemento(nel)
use variaveis
implicit none
integer :: nel,NGP,n1,i,j
real(4) :: E,G,DET,QUOT,SUM
real(4) :: DEE(4,4),SAMP(4,2),COORD(4,2),JAC(2,2),JAC1(2,2)
real(4) :: DER(2,4),DERIV(2,4),BEE(4,8),DBEE(4,8)
real(4) :: BTDB(8,8),KM(8,8),ELD(8),EPS(3),SIGMA(4)
real(4) :: BT(8,4),FUN(4)
integer, parameter :: IDEE=4,IBEE=4,IDBEE=4,IH=4
integer, parameter :: IDOF=8,IBTDB=8,IBT=8,IKM=8
integer, parameter :: IJAC=2,IJAC1=2,NODOF=2,IT=2,IDER=2,IDERIV=2
integer, parameter :: ICOORD=4,NOD=4
integer, parameter :: ISAMP=4

!     INITIALISATION
      E=Ei(imat(nel))
      G=Poisson(imat(nel))
      CALL FMDRAD(DEE,IDEE,E,G)
      NGP=2
      CALL GAUSS(SAMP,ISAMP,NGP)
!     Coordenadas dos nos
      n1=nne*(nel-1)
      COORD(1,1)=x(kon(n1+1))
      COORD(2,1)=x(kon(n1+2))
      COORD(3,1)=x(kon(n1+3))
      COORD(4,1)=x(kon(n1+4))
      COORD(1,2)=y(kon(n1+1))
      COORD(2,2)=y(kon(n1+2))
      COORD(3,2)=y(kon(n1+3))
      COORD(4,2)=y(kon(n1+4))
!!    ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY
      CALL NULL(KM,IKM,IDOF,IDOF)
      do I=1,NGP
        do J=1,NGP
          CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)
          CALL XMATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)
          CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
          CALL XMATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)
          CALL NULL(BEE,IBEE,IH,IDOF)
          CALL FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
          CALL XMATMUL(DEE,IDEE,BEE,IBEE,DBEE,IDBEE,IH,IH,IDOF)
          CALL MATRAN(BT,IBT,BEE,IBEE,IH,IDOF)
          CALL XMATMUL(BT,IBT,DBEE,IDBEE,BTDB,IBTDB,IDOF,IH,IDOF)
          QUOT=SUM*DET*SAMP(I,2)*SAMP(J,2)
          CALL MSMULT(BTDB,IBTDB,QUOT,IDOF,IDOF)
          CALL MATADD(KM,IKM,BTDB,IBTDB,IDOF,IDOF)
        enddo
      enddo

      do i=1,8
        do j=1,8
          RigElem(i,j)=KM(i,j)
        enddo
      enddo

   return
end subroutine Monta_Matriz_Elemento



subroutine Calcula_Incognitas_Secundarias()
use variaveis
implicit none
integer :: NGP,nel,n1,M,k1,i,j,igauss
real(4) :: E,anu,DET,SUM
real(4) :: DEE(4,4),SAMP(4,2),COORD(4,2),JAC(2,2),JAC1(2,2)
real(4) :: DER(2,4),DERIV(2,4),BEE(4,8),DBEE(4,8)
real(4) :: BTDB(8,8),KM(8,8),ELD(8),EPS(3),SIGMA(4)
real(4) :: BT(8,4),FUN(4)
integer, parameter :: IDEE=4,IBEE=4,IDBEE=4,IH=4
integer, parameter :: IDOF=8,IBTDB=8,IBT=8,IKM=8
integer, parameter :: IJAC=2,IJAC1=2,NODOF=2,IT=2,IDER=2,IDERIV=2
integer, parameter :: ICOORD=4,NOD=4
integer, parameter :: ISAMP=4

!     INITIALISATION
      allocate(stress(nelem,4,4))
      NGP=2
      CALL GAUSS(SAMP,ISAMP,NGP)

      do nel=1,nelem
        E=Ei(imat(nel))
        anu=Poisson(imat(nel))
        CALL FMDRAD(DEE,IDEE,E,anu)
!     Coordenadas dos nos
        n1=nne*(nel-1)
        COORD(1,1)=x(kon(n1+1))
        COORD(2,1)=x(kon(n1+2))
        COORD(3,1)=x(kon(n1+3))
        COORD(4,1)=x(kon(n1+4))
        COORD(1,2)=y(kon(n1+1))
        COORD(2,2)=y(kon(n1+2))
        COORD(3,2)=y(kon(n1+3))
        COORD(4,2)=y(kon(n1+4))
        do M=1,4
          k1=ndf*(kon(n1+M)-1)+1
          ELD(2*(M-1)+1)=Desloc(k1)
          ELD(2*(M-1)+2)=Desloc(k1+1)
          !write(io,*)ELD(2*(M-1)+1),'  ',ELD(2*(M-1)+2)
        enddo
!      igauss=0
        do I=1,NGP
          do J=1,NGP
!           igauss=igauss+1
            CALL FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)
            CALL XMATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,IT,NOD,IT)
            CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
            CALL XMATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,IT,IT,NOD)
            CALL NULL(BEE,IBEE,IH,IDOF)
            CALL FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
            CALL MVMULT(BEE,IBEE,ELD,IH,IDOF,EPS)
            CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)
            if (J.eq.1.and.I.eq.1) igauss=2
            if (J.eq.2.and.I.eq.1) igauss=3
            if (J.eq.1.and.I.eq.2) igauss=1
            if (J.eq.2.and.I.eq.2) igauss=4
            stress(nel,igauss,1)=SIGMA(1)
            stress(nel,igauss,2)=SIGMA(2)
            stress(nel,igauss,3)=SIGMA(3)
            stress(nel,igauss,4)=SIGMA(4)
          !write(io,*)nel,' J=',J,' I=',I,SAMP(J,1),SAMP(I,1),stress(nel,igauss,1),stress(nel,igauss,2),stress(nel,igauss,3)
          enddo
        enddo
      enddo

   return
end subroutine Calcula_Incognitas_Secundarias


end module Matriz_Elemento


module Matriz_Global

contains

subroutine Monta_Matriz_Global()
use variaveis
use Matriz_Elemento
implicit none
integer :: i,j,nel

    allocate(RigGlobal(neq,neq))
    do i=1,neq
       do j=1,neq
          RigGlobal(i,j)=0.0  !! Zera matriz global
       enddo
    enddo

    do nel=1,nelem
       call Monta_Matriz_Elemento(nel)
       call Armazena_Elemento_na_Global(nel)
    enddo

   return
end subroutine Monta_Matriz_Global



subroutine Armazena_Elemento_na_Global(nel)
use variaveis
implicit none
integer :: nel
integer :: i,j,i1,j1,i2,j2,n1,n2,k,L,jr,kr,jc,kc

   do i=1,nne
      n1=kon(nne*(nel-1)+i)
      i1=ndf*(i-1)
      j1=ndf*(n1-1)
      do j=1,nne
         n2=kon(nne*(nel-1)+j)
         i2=ndf*(j-1)
         j2=ndf*(n2-1)
         do k=1,ndf
            jr=i1+k
            kr=j1+k
            do L=1,ndf
               jc=i2+L
               kc=j2+L
               RigGlobal(kr,kc)=RigGlobal(kr,kc)+RigElem(jr,jc)
            enddo
         enddo
      enddo
   enddo

   return
end subroutine Armazena_Elemento_na_Global



subroutine Impoe_Condicoes_de_Contorno()
use variaveis
implicit none
integer :: L,no,k1,i,kr,j

   do L=1,nbn
      no=istatus((ndf+1)*(L-1)+1)
      k1=ndf*(no-1)
      do i=1,ndf
         if (istatus((ndf+1)*(L-1)+1+i).eq.1) then
            kr=k1+i
            do j=1,neq
               Carga(j)=Carga(j)-RigGlobal(kr,j)*Prescrito(kr)
               RigGlobal(kr,j)=0.0
               RigGlobal(j,kr)=0.0
            enddo
            RigGlobal(kr,kr)=1.0
            Carga(kr)=Prescrito(kr)
         endif
      enddo
   enddo

   return
end subroutine Impoe_Condicoes_de_Contorno



subroutine Sistema_Linear_Gauss(n,A,B,X,io)
implicit none
integer :: n,io
real(4) :: A(n,n),B(n),X(n)
integer :: i,j,k
real(4) :: c

   do k=1,(n-1)
      c=A(k,k)
      if (abs(c).lt.1.0E-7) then
         write(io,*)'Singularidade na linha ',k
         return
      endif
      do j=1,n
         A(k,j)=A(k,j)/c
      enddo
      B(k)=B(k)/c
      do i=(k+1),n
         c=A(i,k)
         do j=1,n
            A(i,j)=A(i,j)-c*A(k,j)
         enddo
         B(i)=B(i)-c*B(k)
      enddo
   enddo
   if (abs(A(n,n)).lt.1.0E-7) then
      write(io,*)'Singularidade na linha ',k
      return
   endif
   B(n)=B(n)/A(n,n)
   do i=(n-1),1,-1
      do j=(i+1),n
         B(i)=B(i)-A(i,j)*B(j)
      enddo
   enddo
   do i=1,n
      X(i)=B(i)
   enddo

   return
end subroutine Sistema_Linear_Gauss


end module Matriz_Global



program isop4_axis95
use variaveis
use Entrada_Saida
use Matriz_Elemento
use Matriz_Global
implicit none

   call Abre_Arquivos()
   call Entrada_Dados()
   call Monta_Matriz_Global()
   call Impoe_Condicoes_de_Contorno()
   call Sistema_Linear_Gauss(neq,RigGlobal,Carga,Desloc,io)
   call Calcula_Incognitas_Secundarias()
   call Imprime_Resultados()
   call Resultados_GiD()
   call Fecha_Arquivos()

end program isop4_axis95

