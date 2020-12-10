program creat_basinID
! ==================================================
! creat IDs for all basins the number is depend on the catchment area
! of each basin 
! eg: 1 - Amazon, 2 - Congo etc
! 1    Amazon
! 2    Congo
! 3    Mississippi
! 4    Nile
! 5    Parana
! 6    Yenisei
! 7    Ob
! 8    Lena
! 9    Niger
! 10   Amur
! Menaka @ IIS 2019/10/15
! edited to read any river network map
! ==================================================
      implicit none
!
      character*128              ::  fname,buf,camadir,mapname,outdir
      real                       ::  west, east, north, south
      integer                    ::  nx, ny, nflp
      real                       ::  gsize
!
      integer                    ::  ix, iy, k, jx, jy, x, y
      integer                    ::  ud, countnum
! allocatables
      integer,allocatable        ::  nextx0(:,:),  nexty0(:,:) !!  global maps
      real,allocatable           ::  uparea0(:,:)              !!  upstream area 
      !integer,allocatable        ::  nextx(:,:),  nexty(:,:)   !!  regional maps
      integer,allocatable        ::  rivmthX(:), rivmthY(:)    !!  river mouth locations
      integer,allocatable        ::  rivnum(:,:)
      real,allocatable           ::  uparea(:), upareanew(:)!! 
      integer,allocatable        ::  rnk2org(:), org2rnk(:) 
!===========================================================
call getarg(1,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(2,buf)
read(buf,"(A)") mapname
write(*,*) mapname

call getarg(3,buf)
read(buf,"(A)") outdir

!==
fname=trim(camadir)//"/map/"//trim(mapname)//"/params.txt"
print *, fname
open(11,file=fname,form='formatted')
read(11,*) nx
read(11,*) ny
read(11,*) nflp
read(11,*) gsize
read(11,*) west
read(11,*) east
read(11,*) south
read(11,*) north
close(11)
!----
allocate(nextx0(nx,ny),nexty0(nx,ny),uparea0(nx,ny))
allocate(rivmthX(nx*ny),rivmthY(nx*ny))
allocate(rivnum(nx,ny))
!-----
print *, 'read maps'
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nextxy.bin"
print *, fname
open(11,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
read(11,rec=1) nextx0
read(11,rec=2) nexty0
close(11)
!--
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/uparea.bin"
print *, fname
open(11,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
read(11,rec=1) uparea0
close(11)
!-find all river mouth 
countnum=1
do iy=1, ny
    do ix=1, nx
        if (nextx0(ix,iy) == -9 .or. nextx0(ix,iy)== -10) then 
            rivmthX(countnum) = ix
            rivmthY(countnum) = iy
            countnum=countnum+1
        end if 
    end do
end do
!------
countnum=countnum-1
write(*,*) "rivmth:" ,countnum
allocate(uparea(countnum), rnk2org(countnum), org2rnk(countnum), upareanew(countnum))
rivnum=-9999
!--
do k=1,countnum
    x=rivmthX(k)
    y=rivmthY(k)
    !--
    uparea(k)=uparea0(x,y)
end do
!---
call sort_decord(countnum,uparea,upareanew,org2rnk,rnk2org)
write(*,*)org2rnk(1:10)
write(*,*)rnk2org(1:10)
!---
do k=1,countnum
    x=rivmthX(rnk2org(k))
    y=rivmthY(rnk2org(k))
    write(*,*)k,x,y !-180.0 +x * 0.25 , 90.0 - y * 0.25
    rivnum(x,y)=k
    do jx=1, nx
        do jy=1, ny
            if ( nextx0(jx,jy) == -9999 )then !! ocean
                cycle
            end if
            if ( nextx0(jx,jy) == -9 .or. nextx0(jx,jy) == -10)then !! if river mouth or inland termination
                cycle
            end if
            call uord(jx,jy,x,y,nx,ny,nextx0,nexty0,ud) !! if the downstream in the same river
            if (ud == -1) then
              !write(*,*) ud, -180.0 +jx* 0.25 , 90.0 - jy * 0.25
              rivnum(jx,jy)=k
            end if
        end do
    end do
end do
print *, 'write basin ID',trim(mapname)
fname=trim(adjustl(outdir))//"/rivnum_"//trim(mapname)//".bin"
print *, fname
open(21,file=fname,form='unformatted',access='direct',recl=4*nx*ny)
write(21,rec=1) rivnum
close(21)
!--
21 format(i4.4,2x,i4.4,2x,f8.4,2x,f8.4,2x,f8.4)
22 format(a7,2x,a4,2x,a4,2x,a11)
23 format(i7,2x,i4.4,2x,i4.4,2x,e11.4)

fname=trim(adjustl(outdir))//"/rivmth_"//trim(mapname)//".txt"
open(72,file=fname,status='replace')
write(72,22)"Id","lon","lat","uparea"
do k=1, countnum
  write(72,23) k, rivmthX(rnk2org(k)), rivmthY(rnk2org(k)), upareanew(k)
  !write(*,*) k, rivmthX(rnk2org(k)), rivmthY(rnk2org(k)), upareanew(k)
end do
close(72)
deallocate(nextx0,nexty0,uparea0)
deallocate(rivmthX,rivmthY,rivnum)
deallocate(uparea, rnk2org, org2rnk, upareanew)
! ==================================================
end program creat_basinID
!*****************************************************************
subroutine uord(i,j,x,y,nx,ny,nextX,nextY,ud)
implicit none 
!--
! find the (i,j) upstream/downstream of (x,y) 
! upstream ud = -1  downstream ud = +1
! ud = -9999 : another river 
integer                     :: i,j,x,y,nx,ny
integer,dimension(nx,ny)    :: nextX,nextY
!--
integer                     :: ix,iy,iix,iiy,tx,ty,ud
!--
tx=x
ty=y
ix=i
iy=j
ud=-1
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  !write(*,*) iix, iiy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  !write(*,*) "--------------",ix,iy,nextX(ix,iy),nextY(ix,iy),x,y,(nextX(ix,iy)/=tx)! .and. nextY(ix,iy)/=ty)
  if (ix==-9 .or. iy==-9) then ! river mouth
    !write(*,*)"end of the river"
    ud=+1 
    exit
  end if
  if  (ix ==-10 .or. iy ==-10) then ! inland termination 
    !write(*,*)"end of the river"
    ud=+1
    exit
  end if
  if  (ix==-9999 .or. iy==-9999) then
    !write(*,*)"end of the river"
    ud=-9999
    exit
  end if
end do
!---
if (ud==+1) then
  !-
  tx=i
  ty=j
  ix=x
  iy=y
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix ==-9 .or. iy==-9) then ! river mouth
      ud=-9999
      exit
    end if
    if  (ix==-10 .or. iy==-10) then ! inland termination
      !write(*,*)"end of the river"
      ud=-9999
      exit
    end if
    if  (ix==-9999 .or. iy==-9999) then
      !write(*,*)"end of the river"
      ud=-9999
      exit
    end if 
  end do
end if
!write(*,*) ud
return
!---
end subroutine uord 
!*****************************************************************
subroutine sort_decord(n0rec,r1dat,r1out,i1org2rnk,i1rnk2org)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cto   sort an array in decreasing order
!cby   2010/03/31, hanasaki, NIES: H08 ver1.0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
!c parameter (array)
      integer           n0rec
!c parameter (default)
      real              p0mis
      real              p0maxini
      parameter        (p0mis=1.0E20) 
      parameter        (p0maxini=-9.99E20) 
!c in
      real              r1dat(n0rec)     !! original data
!c out
      real              r1out(n0rec)
      integer           i1org2rnk(n0rec) !! original order --> rank
      integer           i1rnk2org(n0rec) !! rank --> original order
!c local
      integer           i0rec
      integer           i0rnk
      integer           i1flg(n0rec)
      real              r0max
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
!c Initialize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      i1flg=0
      r0max=p0maxini
      r1out=p0mis
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Calculate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i0rnk=1,n0rec
        do i0rec=1,n0rec
          if(i1flg(i0rec).eq.0)then
            r0max=max(r1dat(i0rec),r0max)
          end if
        end do
        do i0rec=1,n0rec
          if(r1dat(i0rec).eq.r0max.and.i1flg(i0rec).eq.0)then
            r1out(i0rnk)=r1dat(i0rec)
            i1org2rnk(i0rec)=i0rnk
            i1rnk2org(i0rnk)=i0rec
            i1flg(i0rec)=1
            goto 55
          end if
        end do
 55     continue
        r0max=p0maxini
      end do
      return 
!c
end subroutine sort_decord
!************************************************************************
