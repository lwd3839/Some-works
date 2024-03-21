program DGTD_POST
    
    implicit none
    
    integer :: ii,jj,kk,mm,nn,io,rflag
    doubleprecision :: c0,pi,k0
    complex*16 :: j0,ek
    !-------输入文件_near里读的-------
    integer :: mode_flag,monitor_num
    integer :: num_ext
    doubleprecision :: freq
    doubleprecision,allocatable :: n_ext(:,:)
    complex*16 :: einc,hinc
    complex*16,allocatable :: Jx(:),Jy(:),Jz(:),Jmx(:),Jmy(:),Jmz(:)
    character(len=100) :: filename_char,monitor_id
    
    !------输入文件farpost里读的--------
    doubleprecision :: st_theta,ed_theta,st_phi,ed_phi,dtheta,dphi
    doubleprecision :: itheta,iphi
    integer :: theta_num,phi_num
    doubleprecision :: dir(3)   
    
    doubleprecision :: RCS,ZZ
    complex*16 :: fx,fy,fz,fmx,fmy,fmz,Atheta,Aphi,Ftheta,Fphi
    complex*16 :: etheta,ephi
    
    c0 = 299792458D0
    pi = 3.14159265359
    j0 = (0.0,1.0)
    ZZ = 120*pi
    rflag=0
10010 write(*,*)'Input the filename of farfield (*_near.txt)'
    open(4321,file='post.set')
    read(4321,*)
    read(4321,*)filename_char
    open(1234,file=trim(adjustl(filename_char))//'_near.txt',status='old',iostat=io)
    if(io/=0)then
        write(*,*)'[#] Error occurred while openning file !'
        !write(*,*)' Retry ? [Yes=1]'
        !read(*,*)rflag
        !if(rflag==1)goto 10010
        close(4321)
        stop
    endif
    
    write(*,*)'Input the angle theta of post (start,end,number) [degree]'
    read(4321,*)
    read(4321,*)st_theta,ed_theta,theta_num
    dtheta = (ed_theta-st_theta)/theta_num
    if(theta_num>1)then
        dtheta = (ed_theta-st_theta)/(theta_num-1)
    else
        dtheta = (ed_theta-st_theta)/theta_num
    endif
    write(*,*)'Input the angle phi of post (start,end,number) [degree]'
    read(4321,*)
    read(4321,*)st_phi,ed_phi,phi_num
    close(4321)
    if(phi_num>1)then
        dphi = (ed_phi-st_phi)/(phi_num-1)
    else
        dphi = (ed_phi-st_phi)/phi_num
    endif
    !角度的定义与feko相同
    
    read(1234,*) mode_flag
    read(1234,*) monitor_num
    do ii = 1,monitor_num
        read(1234,*)num_ext
        read(1234,*)freq
        read(1234,*)einc
        read(1234,*)hinc
        allocate(n_ext(3,num_ext))
        allocate(Jx(num_ext),Jy(num_ext),Jz(num_ext))
        allocate(Jmx(num_ext),Jmy(num_ext),Jmz(num_ext))
        do jj = 1,num_ext
            read(1234,*)n_ext(1,jj)
            read(1234,*)n_ext(2,jj)
            read(1234,*)n_ext(3,jj)
            read(1234,*)Jx(jj)
            read(1234,*)Jy(jj)
            read(1234,*)Jz(jj)
            read(1234,*)Jmx(jj)
            read(1234,*)Jmy(jj)
            read(1234,*)Jmz(jj)
        enddo
        write(*,*)' File reading finished!'
        
        write(monitor_id,*)ii
        open(2345,file=trim(adjustl(filename_char))//'_monitor_'//trim(adjustl(monitor_id))//'_far.txt')
        write(2345,*)' theta[deg] phi[deg] rEtheta iEtheta rEphi iEphi RCS(sm) RCS(dBsm)'
        
        k0 = 2*pi*freq/c0

        itheta = st_theta-dtheta
        iphi = st_phi-dphi
        do jj = 1,theta_num
            itheta = itheta + dtheta
            do kk = 1,phi_num
                iphi = iphi + dphi
                dir(1) = sind(itheta)*cosd(iphi)    !x
                dir(2) = sind(itheta)*sind(iphi)    !y
                dir(3) = cosd(itheta)               !z
                fx = (0,0);fmx = (0,0)
                fy = (0,0);fmy = (0,0)
                fz = (0,0);fmz = (0,0)                           

                do mm = 1,num_ext
                    !ek = exp(j0*dot_product(dir,n_ext(1:3,mm)))
                    !ek = exp(j0*dot_product(dir,n_ext(1:3,mm)))!/4/pi
                    ek = exp(j0*k0*dot_product(dir,n_ext(1:3,mm)))
                    fx = fx+Jx(mm)*ek
                    fy = fy+Jy(mm)*ek
                    fz = fz+Jz(mm)*ek
                    fmx = fmx+Jmx(mm)*ek
                    fmy = fmy+Jmy(mm)*ek
                    fmz = fmz+Jmz(mm)*ek
                enddo
                
                !A
                Atheta = fx*cosd(itheta)*cosd(iphi)+fy*cosd(itheta)*sind(iphi)-fz*sind(itheta)
                Aphi = -1*fx*sind(iphi)+fy*cosd(iphi)
                !F
                Ftheta = fmx*cosd(itheta)*cosd(iphi)+fmy*cosd(itheta)*sind(iphi)-fmz*sind(itheta)
                Fphi = -1*fmx*sind(iphi)+fmy*cosd(iphi)
                
                Etheta = -1*j0*k0*exp(-1*j0*k0)/(4*pi)*(120*pi*Atheta+Fphi)
                Ephi = j0*k0*exp(-1*j0*k0)/(4*pi)*(-1*120*pi*Aphi+Ftheta)
                
                !Etheta = j0*k0/4/PI*((fmx*SIND(iphi)-fmy*cosD(iphi))-zz*(fx*cosD(itheta)*cosD(iphi)+ &
                !    fy*cosD(itheta)*sinD(iphi)-fz*sIND(itheta)))
                !Ephi =  j0*k0/4/PI*(zz*(fx*sIND (iphi)-fy*CosD(iphi))+(fmx*cosD(itheta)*cosD(iphi)+ &
                !    fmy*cosD(itheta)*SIND(iphi)-fmz*sIND(itheta)))
                
                !Etheta = -1*exp(-1*j0*k0)*j0*k0*(120*pi*Atheta+Fphi)
                !Ephi = -1*exp(-1*j0*k0)*j0*k0*(120*pi*Aphi-Ftheta)
                
                RCS = 10*log10(4*pi*(abs(Etheta)**2+abs(Ephi)**2)/cdabs(Einc)**2)
                
                write(2345,'(8(E15.6,X))') itheta,iphi,real(Etheta),imag(Etheta),real(Ephi),imag(Ephi),4*pi*(abs(Etheta)**2+abs(Ephi)**2),RCS
            enddo
        enddo
        close(2345)
        deallocate(Jx,Jy,Jz,Jmx,Jmy,Jmz,n_ext)
    enddo
    close(1234)
    
endprogram
