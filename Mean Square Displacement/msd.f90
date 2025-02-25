program msd
! initialising variables
character(len=10) :: a,b
integer :: c, num_frames = 25000 , num_mol = 884
real , allocatable :: r(:,:,:) , box_len(:), msd(:)
real, allocatable :: r1(:,:) ! dummy 
integer :: dt
double precision,allocatable :: sd(:,:),msd1(:)
double precision :: t1, t2 

! allocating
allocate(r(num_mol,num_frames,3) , box_len(num_frames))
allocate(msd(num_frames),r1(num_frames,3))
allocate(sd(num_frames,3),msd1(num_frames))

print*, " You have", num_frames, " number of frames"
print*, " You have", num_mol , " number of water molecules"
call cpu_time(t1)   ! print value in secondds

!opening of input file
open(unit = 10, file = "./all_oxygen.gro", status = "old", action ="read") ! input trajectory file

! reading input files i.e. traj file of O
r(:,:,:) = 0
box_len(:) = 0
do i = 1 , num_frames
        read(10,*)
        read(10,*)
        do j = 1, num_mol
                read(10,*) a,b,c, r(j,i,1), r(j,i,2), r(j,i,3)
                write(100,*) a,b,c,r(j,i,1), r(j,i,2), r(j,i,3)
        end do
        read(10,*) box_len(i)
        write(100,*) box_len(i)
end do

print*, "Input trajectory file read sucessfully, Now unfolding PBC applied"

! Unfolding PBC
do j = 1, num_mol
        do i = 1, num_frames-1 ! -1 because at last step it will go out of bounds
                do k = 1,3
                        r(j,i+1,k) = r(j,i+1,k) - r(j,i,k)
                        r(j,i+1,k) = r(j,i+1,k) - (anint(r(j,i+1,k)/box_len(i)) * box_len(i))
                        r(j,i+1,k) = r(j,i+1,k) + r(j,i,k)
                end do
        end do
end do

! checking removing of pbc
do j = 1, num_mol
        do i = 1, num_frames
                write(400,*) i , r(j,i,:)
        end do
end do


! Now calculating MSD
print*, " Calculating MSD, sit back and relax"
msd(:)  = 0.0
msd1(:) = 0.0
r1(:,:) = 0.0
sd(:,:) = 0.0
do i = 1, num_mol
        r1 = r(i,:,:)
outer : do l = 1, num_frames
!        msd1(l) = 0.0
!        sd(l,:)  = 0.0
        inner : do k = 1, num_frames
                if ( k+l .gt. (num_frames)) exit inner ! without it, diffrence will go upto (ntime+2)-(ntime+1)
                sd(l,:)  = sd(l,:) + ((r1(k+l,:) - r1(k,:))**2)
                end do inner
                msd1(l) = (sum(sd(l,:)))/real(num_frames-l)
        end do outer
        if (mod(i,20) ==0) then
                print*, i
        end if

end do
msd = msd1 / dble(num_mol)

deallocate(r,r1)

print*, "MSD is calculated, output file is fort.500"

do dt = 1, num_frames-1
        write(500,*) dt , msd(dt)
end do

print*, "Everything is okay, Well done"
call cpu_time(t2)
write(*,'(A , F6.3)') " Total time taken for code running is (in minutes)" , (t2-t1)/60


end program msd
