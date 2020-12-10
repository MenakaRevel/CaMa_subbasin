program sub_basin
!$ use omp_lib
! Find the convolvation points => if two braches carrying more than 40% of upstream area
! allocate river pixels flowing through the covolvation points
! 10% of total catchment area of the river
! minmum area as input threshold
! Menaka@IIS 2019/10/22
implicit none
character*128                   :: fname,buf,camadir,mapname,outdir,header
!character*8                     :: yyyymmdd,nxtyyyymmdd
!real                            :: assimN,assimS,assimW,assimE
real                            :: lat,lon
!character*2                     :: swot_day
!real,allocatable                :: swot_obs(:,:),global_xa(:,:,:)
integer*4                       :: ix,iy,jx,jy,i,kx,ky !lon_cent,lat_cent,patch_size,patch_side,,countnum,patch_nums,j
!integer,parameter               :: latpx=720,lonpx=1440
integer                         :: nx, ny, nflp
real                            :: gsize
real                            :: west, east, north, south
real,allocatable                :: rivwth(:,:),nextdst(:,:),rivlen(:,:),uparea(:,:),subbasin(:,:),subriver(:,:),basin(:,:) !,gauss_weight,
integer,allocatable             :: nextX(:,:),nextY(:,:),ocean(:,:),rivseq(:,:),rivnum(:,:),targetp(:,:)

!real,allocatable                :: Wvec(:),lag(:),local_lag(:)!global_sum_xam(:,:),
!real,allocatable                :: wgt(:),local_wgt(:)
!real                            :: wt,lag_dist,conflag!lag_distance,Gauss_wt,
!integer*4                       :: i_m,j_m
integer,allocatable             :: rivid(:),xlist(:),ylist(:)
real,allocatable                :: alist(:)
!integer*4                       :: target_pixel,fn
!character*8                     :: llon,llat
real                            :: perct,threshold,catchment_area,area_min,upbarea
integer                         :: ios,countp,ud!,rivernum
integer,allocatable             :: xconlist(:),yconlist(:),xconl(:),yconl(:)
real,allocatable                :: convarea(:),areal(:),arealnew(:)!! for decreasing order
integer,allocatable             :: rnk2org(:),org2rnk(:) 
real                            :: numm,lag_dist, jbsn
integer                         :: N,num,seq
integer,dimension(8)            :: xl,yl
real,dimension(8)               :: al
!--
write(*,*) "create subbasin"
print*, "######################################################"
print*, "create subbasin"
print*, "USAGE: Input Uparea Threshold [km2]"
print*, "./subbasin $Map $CaMadir $Outdir $UPA_THRS $PRCT"
! call getarg(1,buf)
! read(buf,*) N ! the number for paticular river

call getarg(1,buf)
read(buf,"(A)") mapname
write(*,*) mapname

call getarg(2,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(3,buf)
read(buf,"(A)") outdir

call getarg(4,buf)
read(buf,*) threshold

call getarg(5,buf)
read(buf,*) perct

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
!patch_side=patch_size*2+1
!patch_nums=patch_side**2


! 21 format(i4.4,2x,i4.4,2x,f10.7)
! 22 format(a4,2x,a4,2x,a8)
! 23 format(i4.4,2x,i4.4,2x,i4.4,2x,e11.4)
24 format(a10,2x,a4,2x,a4,2x,a11,2x,a11)
25 format(f10.3,2x,i4.4,2x,i4.4,2x,e11.4,2x,e11.4)
!$ write(*,*)"omp threads",omp_get_num_threads()

allocate(rivwth(nx,ny),nextdst(nx,ny),rivlen(nx,ny),uparea(nx,ny),subbasin(nx,ny),subriver(nx,ny),basin(nx,ny))
allocate(nextX(nx,ny),nextY(nx,ny),ocean(nx,ny),rivseq(nx,ny),rivnum(nx,ny),targetp(nx,ny))


! read river width
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivwth.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth"
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nextxy.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextX
    read(34,rec=2) nextY
else
    write(*,*) "no file nextXY at:",fname
end if
close(34)

! make ocean mask from nextX data (1 is ocean; 0 is not ocean)
ocean = (nextX==-9999) * (-1)

! read river length
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)
! read distance to next grid
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
else
    write(*,*) "no file nextdst",fname
end if
close(34)
! read rivseq file
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/rivseq.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivseq
else
    write(*,*) "no file rivseq",fname
end if
close(34)
! read uparea file
fname=trim(adjustl(camadir))//"/map/"//trim(mapname)//"/uparea.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) uparea
else
    write(*,*) "no file uparea",fname
end if
close(34)
! read rivnum file
fname=trim(adjustl(outdir))//"/rivnum_"//trim(mapname)//".bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivnum
else
    write(*,*) "no file rivnum",fname
end if
close(34)
!--get N
N=100 !int(maxval(rivnum))
!--
allocate(xconlist(100000),yconlist(100000),convarea(100000))
!allocate(xlist(11*11),ylist(11*11),area(11*11))
!--
xconlist=-9999
yconlist=-9999
convarea=-9999
targetp=0
! read rivmth 
allocate(rivid(N),xlist(N),ylist(N),alist(N))
fname=trim(adjustl(outdir))//"/rivmth_"//trim(mapname)//".txt"
open(34,file=fname,form="formatted",iostat=ios)
if (ios==0) then
     read(34,*) header
     do i=1,N
         read(34,*) rivid(i),xlist(i),ylist(i),alist(i)
     end do
else
    write(*,*) "no rivmth", fname
end if
close(34)
! write(*,*)header
! write(*,23)rivid(1),xlist(1),ylist(1),alist(1)
! write(*,23)rivid(2),xlist(2),ylist(2),alist(2)
! write(*,23)rivid(3),xlist(3),ylist(3),alist(3)
! write(*,23)rivid(4),xlist(4),ylist(4),alist(4)
!--
do iy=1, ny
  do ix=1, nx
    if( uparea(ix,iy)>0 ) uparea(ix,iy)=uparea(ix,iy)*1.e-6
  end do
end do
!--
subbasin=real(rivnum) !-9999.0 !
subriver=-9999.0 !real(rivnum)*uparea>=
basin=0.0
fname=trim(adjustl(outdir))//"/outsub_"//trim(mapname)//".txt"
open(72,file=fname,form="formatted",iostat=ios)
write(72,24)"Id","lon", "lat", "uparea[km2]", "subarea"
!--
! do parallel
!$omp parallel default(none)&
!$omp& private(lon_cent,lat_cent,lat,lon,llon,llat,fname,fn,lag,countnum,i,j,i_m,j_m,threshold)&
!$omp& shared(ocean,rivwth,targetp,countp,patch_size,weightage,gauss_weight)
!$omp do
!--
do num=1, N !max(rivnum)
    countp=1
    catchment_area=alist(num)*1.e-6
    xconlist(1)=xlist(num)
    yconlist(1)=ylist(num)
    convarea(1)=alist(num)*1.e-6
    !--
    !write(*,*) num
    !print*, "river :", num
    !--
    !write(*,*)catchment_area,threshold
    !area_min=min(0.1*catchment_area,threshold)
    !!!!! for more uniform catchments / for extream cases use total basin area [very small basins]
    area_min=min(catchment_area,threshold)
    !subriver=real(rivnum)*(rivnum==num)*(uparea>=area_min)*-1.0
    !write(*,*) xconlist(1), yconlist(1)
    do ix = 1,nx !
        do iy = 1, ny !
            lat = 90.0 - (iy-1.0)*gsize
            lon = (ix-1.0)*gsize - 180.0
            !remove ocean
            if (ocean(ix,iy)==1) then
                cycle
                !continue
            end if
            ! remove rivwth <= 0m
            if (rivwth(ix,iy)<=0.0) then
                cycle
            end if
          ! remove other rivers
            if (rivnum(ix,iy)/=num) then
                cycle
            end if
            !!!!!! old ****** consider only upstream area >= 10% * total catchment area of river or thershold
            !consider only upstream area >= thershold
            !uparea(ix,iy)=uparea(ix,iy)*1.e-6
            if (uparea(ix,iy)<=area_min) then
                cycle
            end if
            !write(*,*)num, xlist(num),ylist(num)!ix,iy! "===",lat,lon,"==="
            ! river mouth , inland termination point
            if (nextX(ix,iy)==-9 .or. nextX(ix,iy)==-10) then
                !xconlist(countp)=lon
                !yconlist(countp)=lat
                !convarea(countp)=uparea(ix,iy)
                !countp=countp+1
                !write(*,*)ix,iy
                cycle
            end if
            !! remove small sub basins area > area_min
            !if (uparea(ix,iy) < threshold) then
            !    cycle
            !end if
            !--
          !  call find_convolve(ix,iy,area_min,perct,nx,ny,nextX,nextY,uparea,k,xl,yl,al)
          !  !-- consider only two or more tributaries
          !  if (k/=2) then
          !      !write(*,*) "K<2"
          !      cycle
          !      !write(*,*) "K<2"
          !  end if
          !  do i=1,k
          !      xconlist(countp+i)=xl(i)
          !      yconlist(countp+i)=yl(i)
          !      convarea(countp+i)=al(i)
          !      !write(*,*)"*********************"
          !      !write(*,*) i,countp+i,xconlist(countp+i),xl(i)!,yl(i),al(i)
          !  end do
            !###########
            jx=nextX(ix,iy)
            jy=nextY(ix,iy)
            !uparea(jx,jy) = uparea(jx,jy)*1.e-6
            if( uparea(jx,jy)>uparea(ix,iy)+area_min )then
              !print*, "selected :",countp+1,ix, iy
              xconlist(countp+1)=ix
              yconlist(countp+1)=iy
              convarea(countp+1)=uparea(ix,iy)
              countp=countp+1
            end if
            !--
            !write(*,*)"=====================", num,countp,k!,countp+k,xconlist(1)
            !write(*,*)num, xlist(num),ylist(num)!ix,iy! "===",lat,lon,"==="
            !countp=countp+1
        end do
    end do
    !---
    !write(*,*) xconlist(1:5)
    !countp=countp!+10
    !write(*,*) countp!, k
    if (countp <=1) then
        cycle
    end if
    allocate(xconl(countp),yconl(countp),areal(countp),arealnew(countp),org2rnk(countp),rnk2org(countp))
    xconl=-9999
    yconl=-9999
    areal=-9999.0
  !--- 
    xconl=xconlist(1:countp)
    yconl=yconlist(1:countp)
    !areal=convarea(1:countp)
    !arealnew=areal
    !print*, xconl
    print*, "river :", num, countp
    print*, "============================================"
    areal(1)=0.0d0
    do i=2, countp
        !print*, areal(i), xconl(i), xconlist(i)!, yconl(i)
        call lag_distance(xconl(i),yconl(i),xconl(1),yconl(1),nx,ny,nextX,nextY,nextdst,lag_dist)
        areal(i)=lag_dist+mod(real(i),2.0)*0.01
        !write(*,*) i,areal(i),xconl(i),yconl(i)
    end do
    !--
    !do i=1, countp
    !    write(*,*)  areal(i)!,areal(rnk2org(i)), xconl(rnk2org(i))
    !end do
    arealnew=areal
    !---
    !deallocate(xconlist,yconlist,convarea)
    !---
    !call sort_decord(countp,areal,arealnew,org2rnk,rnk2org)
    call sort_incord(countp,areal,arealnew,org2rnk,rnk2org)
    !--
    !do i=1, countp
    !    write(*,*)  areal(i), arealnew(i)!,rnk2org(i),org2rnk(i), xconl(rnk2org(i))
    !end do
    !do i=1, countp
    !    write(*,*)  arealnew(i),areal(rnk2org(i)), xconl(rnk2org(i))
    !end do
    do i=1, countp
      numm=real(num)+real(i-1)*0.001
      subbasin(xconl(rnk2org(i)),yconl(rnk2org(i)))=numm
      basin(xconl(rnk2org(i)),yconl(rnk2org(i)))=numm
      !write(*,*)  arealnew(i),areal(rnk2org(i)), xconl(rnk2org(i))
    end do
    !---
    !num=real(rivernum)
    numm=0.0
    do i=1,countp
        !write(*,*) i, countp,arealnew(i), threshold
        !-- consider only uparea > area_min
        !if (arealnew(i) < threshold) then
        !    cycle
        !end if
        numm=real(num)+real(i-1)*0.001
        !write(*,*) numm, i, countp,arealnew(i), threshold
        !write(72,25)numm,xconl(rnk2org(i)),yconl(rnk2org(i)),uparea(xconl(rnk2org(i)),yconl(rnk2org(i)))
        !--
        do ix = 1,nx !int((assimW+180)*4+1),int((assimE+180)*4+1),1
            do iy = 1, ny !int((90-assimN)*4+1),int((90-assimS)*4+1),1
                lat = 90.0-(iy-1.0)*gsize
                lon = (ix-1.0)*gsize-180.0
                ! remove other rivers
                if (rivnum(ix,iy)/=num) then
                    cycle
                end if
                !remove ocean
                if (ocean(ix,iy)==1) then
                    cycle
                    !continue
                end if
                ! remove rivwth <= 0m
                if (rivwth(ix,iy)<=0.0) then
                    cycle
                end if
                !--
                !numm=real(num)+real(i-1)*0.001
                !-- consider only upstream pixels
                !write(*,*) org2rnk(i)
                !call uord(ix,iy,xconl(org2rnk(i)),yconl(org2rnk(i)),lonpx,latpx,nextX,nextY,ud)
                !  
                if( nextX(ix,iy)>0 .and. basin(ix,iy)==0.0) then
                  jx=ix
                  jy=iy
                  do while( nextX(jx,jy)>0 .and. basin(jx,jy)==0.0 )
                    kx=nextx(jx,jy)
                    ky=nexty(jx,jy)
                    jx=kx
                    jy=ky
                  end do
                  jbsn=basin(jx,jy)
                  print*, jx,jy,jbsn
                  !----------------
                  jx=ix
                  jy=iy
                  do while( nextX(jx,jy)>0 .and. basin(jx,jy)==0.0 )
                    basin(jx,jy)=jbsn
                    subbasin(jx,jy)=jbsn
                    kx=nextx(jx,jy)
                    ky=nexty(jx,jy)
                    jx=kx
                    jy=ky
                  end do 
                end if
            end do
        end do
      end do
      !deallocate(xconlist,yconlist,convarea)
      ! write outsub.txt
      numm=0.0
      do i=1,countp
        upbarea=0.0 
        seq=1000000
        numm=real(num)+real(i-1)*0.001
        do ix=1,nx !int((assimW+180)*4+1),int((assimE+180)*4+1),1
          do iy=1,ny
            if (subbasin(ix,iy)/=numm) cycle
            if (rivseq(ix,iy)==1) cycle
            if (uparea(ix,iy)<threshold) cycle
            if (rivseq(ix,iy)<seq) then
              upbarea=uparea(ix,iy)
            end if
          end do
        end do
        write(72,25)numm,xconl(rnk2org(i)),yconl(rnk2org(i)),uparea(xconl(rnk2org(i)),yconl(rnk2org(i))),uparea(xconl(rnk2org(i)),yconl(rnk2org(i)))-upbarea
      end do
     deallocate(xconl,yconl,areal,arealnew,org2rnk,rnk2org) 
end do
!$omp end do
!$omp end parallel


deallocate(rivid,xlist,ylist,alist)
deallocate(xconlist,yconlist,convarea)
!--sub basin--
fname=trim(adjustl(outdir))//"/subbasin_"//trim(mapname)//".bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) subbasin
else
    write(*,*) "no subbasin", fname
end if
close(84)
!--sub river--
fname=trim(adjustl(outdir))//"/subriver_"//trim(mapname)//".bin"
open(84,file=fname,form="unformatted",access="direct",recl=4*nx*ny,status="replace",iostat=ios)
if(ios==0)then
    write(84,rec=1) subriver
else
    write(*,*) "no subriver", fname
end if
close(84)
close(72)
deallocate(rivwth,nextdst,rivlen,uparea,subbasin,subriver,basin)
deallocate(nextX,nextY,ocean,rivseq,rivnum,targetp)
!allocate(xconl,yconl,areal,arealnew)
!allocate(rivid,xlist,ylist,alist)
end program sub_basin
!*****************************************************************
subroutine find_convolve(x,y,armin,perct,nx,ny,nextX,nextY,uparea,k,xconlist,yconlist,area)
implicit none 
!--
integer                             :: k,i,j,x,y,nx,ny,ix,iy
real                                :: armin,perct ! precentage of upstream to define confulence
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: uparea
!
integer,dimension(8)                :: xconlist,yconlist
real,dimension(8)                   :: area
real,parameter                      :: threshold=1.0d0 ! 30% of area min
!--
xconlist=-9999
yconlist=-9999
area=-9999.0
k=1
do i=x-1,x+1 
    do j=y-1,y+1
        call ixy2iixy(i,j, nx, ny, ix, iy)
        ! next pixle same as x,y
        if (nextX(ix,iy)==x .and. nextY(ix,iy)==y) then  
            if (uparea(ix,iy)>=threshold*armin) then
          ! consider the precntage of contribution
                if (uparea(ix,iy)>=perct*uparea(x,y)) then  
                    !write(*,*) k,uparea(ix,iy)
                    xconlist(k)=ix
                    yconlist(k)=iy
                    area(k)=uparea(ix,iy)
                    !write(*,*) k, xconlist(k),yconlist(k),area(k)
                    k=k+1
                end if
            end if
        end if
    end do
end do
!--
k=k-1
!---
return
!---
end subroutine find_convolve 
!*****************************************************************
! subroutine upbasin_area(x,y,countp,nx,ny,xconl,yconl,uparea,rivseq,area)
! implicit none
! !---
! integer                             :: k,i,j,x,y,nx,ny,ix,iy,countp,ud
! real                                :: armin,perct ! precentage of upstream to define confulence
! integer,dimension(countp)           :: xconl,yconl
! integer,dimension(nx,ny)            :: nextX,nextY,rivseq
! real,dimension(nx,ny)               :: uparea
! !----
! do i=1,countp
!     ix=xconl(i)
!     iy=yconl(i)
!     call uord(ix,iy,x,y,nx,ny,nextX,nextY,ud)
!     if (ud==-1) then

!*****************************************************************
subroutine find_convolve_old(x,y,armin,perct,nx,ny,nextX,nextY,uparea,k,xconlist,yconlist,area)
implicit none 
!--
integer                             :: k,i,j,x,y,nx,ny,ix,iy
real                                :: armin,perct ! precentage of upstream to define confulence
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: uparea
!
integer,dimension(11*11)            :: xconlist,yconlist
real,dimension(11*11)               :: area
real,parameter                      :: threshold=1.0d0 ! 30% of area min
!--
xconlist=-9999
yconlist=-9999
area=-9999.0
k=1
do i=x-5,x+5 
    do j=y-5,y+5
        call ixy2iixy(i,j, nx, ny, ix, iy)
        ! next pixle same as x,y
        if (nextX(ix,iy)==x .and. nextY(ix,iy)==y) then  
            if (uparea(ix,iy)>=threshold*armin) then
                ! consider the precntage of contribution
                if (uparea(ix,iy)>=perct*uparea(x,y)) then  
                     !write(*,*) k,uparea(ix,iy)
                     xconlist(k)=ix
                     yconlist(k)=iy
                     area(k)=uparea(ix,iy)
                     !write(*,*) k, xconlist(k),yconlist(k),area(k)
                     k=k+1
                end if
            end if
        end if
    end do
end do
!--
k=k-1
!---
return
!---
end subroutine find_convolve_old
!*****************************************************************
subroutine uord(i,j,x,y,nx,ny,nextX,nextY,ud)
implicit none 
!--
! find the (i,j) upstream/downstream of (x,y) 
! upstream ud = -1  downstream ud = +1
! ud = -9999 : another river 
integer                     :: i,j,x,y,nx,ny
integer,dimension(nx,ny)    :: nextX,nextY
!real,dimension(nx,ny)       :: rivlen
!--
integer                     :: ix,iy,iix,iiy,tx,ty,ud !pixel,
!--
tx=x
ty=y
ix=i
iy=j
ud=-1
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  !write(*,*) "--------------",ix,iy,nextX(ix,iy),nextY(ix,iy),x,y,(nextX(ix,iy)/=tx)! .and. nextY(ix,iy)/=ty)
  if (ix==-9 .or. iy==-9) then
    ud=+1 
    exit
  end if
  !write(*,*)"end of the river"
  if (ix == -10 .or. iy == -10) then  ! inland termination
    ud=+1 
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
    !write(*,*)iix,iiy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix ==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if 
    if (ix == -10 .or. iy == -10) then ! inland termination
      ud=-9999
      exit
    end if 
  end do
end if
return
!---
end subroutine uord 
!****************************************************************
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
                 !write(*,*)  r1dat(i0rec)
                 goto 55
             end if
         end do
 55      continue
         r0max=p0maxini
     end do
     return 
!c
end subroutine sort_decord
!************************************************************************
subroutine sort_incord(n0rec,r1dat,r1out,i1org2rnk,i1rnk2org)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cto   sort an array in increasing order
!cby   2010/03/31, hanasaki, NIES: H08ver1.0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none
!c parameter (array)
integer                                        :: n0rec
!c parameter (default)
real,parameter                                 :: p0mis=1.0E20
real,parameter                                 :: p0minini=9.99E20
!c in
real,dimension(n0rec)                          :: r1dat     !! original data
!c out
real,dimension(n0rec)                          :: r1out
integer,dimension(n0rec)                       :: i1org2rnk !! original order --> rank
integer,dimension(n0rec)                       :: i1rnk2org !! rank --> original order
!c local
integer                                        :: i0rec            !! index of record
integer                                        :: i0rnk            !! index of rank
integer,dimension(n0rec)                       :: i1flg(n0rec)
real                                           :: r0min
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Initialize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
i1flg=0
r0min=p0minini
r1out=p0mis
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c Calculate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
     do i0rnk=1,n0rec 
       do i0rec=1,n0rec
         !write(*,*) r1dat(i0rec)
         if(i1flg(i0rec).eq.0)then
           r0min=min(r1dat(i0rec),r0min)
         end if
       end do
       do i0rec=1,n0rec
         if(r1dat(i0rec).eq.r0min.and.i1flg(i0rec).eq.0)then
           !write(*,*) r1dat(i0rec) 
           r1out(i0rnk)=r1dat(i0rec)
           i1org2rnk(i0rec)=i0rnk
           i1rnk2org(i0rnk)=i0rec
           i1flg(i0rec)=1
           !write(*,*) r1out(i0rnk),r1dat(i0rec),i0rnk, i0rec
           goto 55
         end if
       end do
 55    continue
       r0min=p0minini
     end do
     return
!c
end subroutine sort_incord
!*****************************************************************
subroutine lag_distance(i,j,x,y,nx,ny,nextX,nextY,nextdst,lag_dist)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
real                                :: lag_dist
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length,rl
!--
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
!--
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  lag_dist=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-10 .or. iy==-10) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!--
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    !---
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-10 .or. iy==-10) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !-- half of the present grid
    rl=anint((nextdst(ix,iy)/1000.0)*100)/100.0
    length=length+rl!/2.0
  end do
end if
!-- 
if (ud==-9999) then
  lag_dist=-9999
elseif (ud==0) then
  lag_dist=0.0
else
  lag_dist=length
end if
!---
return
!---
end subroutine lag_distance
!*****************************************************************
subroutine wgt_consistancy(i,j,x,y,nx,ny,nextX,nextY,weightage,thersold,conflag)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: weightage
!--
real                                :: conflag,thersold
integer                             :: ix,iy,iix,iiy,tx,ty,ud
real                                :: length!,rl
!--
if (i==x .and. j==y) then
  ud=0
else
  ud=-1
end if
!--
if (ud==-1) then
  tx=x
  ty=y
  ix=i
  iy=j
  length=0.0
  conflag=1.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy 
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix==-9 .or. iy==-9) then
      ud=+1
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=+1
      exit
    end if
    !--compare the weightage and thersold
    if (weightage(ix,iy) < thersold) then
      conflag=0.0
    end if 
  end do
end if
!--
if (ud==+1) then
  tx=i
  ty=j
  ix=x
  iy=y
  length=0.0
  conflag=1.0
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    !---
    if (ix==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if
    if (ix==-9999 .or. iy==-9999) then
      ud=-9999
      exit
    end if
    !--compare the weightage and thersold
    if (weightage(ix,iy) < thersold) then
      conflag=0.0
    end if 
  end do
end if
!-- 
if (ud==-9999) then
  conflag=0.0
elseif (ud==0) then
  conflag=1.0
else
  conflag=conflag
end if
!---
return
!---
end subroutine wgt_consistancy
!**************************************************
subroutine read_wgt(fname,nx,ny,weightage)
!$ use omp_lib    
implicit none
character*128                      :: fname 
integer                            :: fn,nx,ny,ios
real,dimension(nx,ny)              :: weightage
fn=34
!$ fn= fn + omp_get_thread_num()
!!$ write(*,*) fn
open(fn,file=fname,form="unformatted",access="direct",recl=4*ny*nx,status="old",iostat=ios)
if(ios==0)then
    read(fn,rec=1) weightage
else
    write(*,*) "no weightage", fname
end if
close(fn)
!--
return
!---
end subroutine read_wgt
!**************************************************
function Gauss_wt(lag)
implicit none
real                                :: lag,Gauss_wt
real,parameter                      :: sigma=1000.0 !1000 km 
!---
Gauss_wt=exp(-(lag**2.0/(2.0*sigma**2.0)))  
!---
return
!---
end function Gauss_wt
!***************************************************   
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     ix, nx
!-- for output ----------
integer                     roundx
!------------------------
if (ix .ge. 1) then
  roundx = ix - int((ix -1)/nx)*nx
else
  roundx = nx - abs(mod(ix,nx))
end if 
return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
implicit none
!- for input -----------------
integer                   ix, iy, nx, ny
!- for output ----------------
integer                   iix, iiy,roundx
!-----------------------------
if (iy .lt. 1) then
  iiy = 2 - iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else if (iy .gt. ny) then
  iiy = 2*ny -iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else
  iiy = iy
  iix = roundx(ix, nx)
end if
return
end subroutine ixy2iixy
!*****************************************************************        
