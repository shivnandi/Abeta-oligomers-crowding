! This code will calculate tetrahedral order  parameter of water molecule
program qtet

character(len=50) :: a
integer :: num_frames , num_mol 
real , allocatable :: r(:,:,:) , box_len(:)
real, allocatable :: r1(:,:)! dummy
real, allocatable :: dist(:,:), r2(:,:,:) ! another dummy
real :: temp_r2(3)
real :: sum_1, sum_qtet, total_qtet, theta, sum2
double precision :: t1, t2
real, allocatable :: qtet1(:,:), qtet2(:)
character(len=50) :: file1

! input from user
!write(*,*) 'Enter number of frames'
!read(*,*) num_frames
!
!write(*,*) 'Enter number of water molecules'
!read(*,*) num_mol
!
!write(*,*) 'Enter traj gro file'
!read(*,*) file1

num_frames = 1
num_mol = XXX

! opening data file
open(unit=10, file = 'heavy.gro' , status="old", action="read") ! input file


! allocating
allocate(r(num_frames,num_mol,3) , box_len(num_frames))
allocate(r1(num_mol,3))
allocate(dist(num_mol,num_mol), r2(num_mol,num_mol,3))
allocate(qtet1(num_frames,num_mol), qtet2(num_frames))
! reading input files i.e. traj file of O
print*, "Thanks, let me handle rest"
r(:,:,:) = 0
box_len(:) = 0
do i = 1 , num_frames
        read(10,*)
        read(10,*)
        do j = 1, num_mol
                read(10,'(1a23,3f7.3)') a, r(i,j,1), r(i,j,2), r(i,j,3)
               ! write(100,'(1a23,3f7.3)') a,r(i,j,1), r(i,j,2), r(i,j,3)
        end do
        read(10,*) box_len(i)
        !write(100,*) box_len(i)
end do
close(10)

total_qtet = 0.0

! MIC and calculating disatnce between all the atoms and saving it in a variable
print*, "Input file read successfully, Sit back and relax"
qtet1(:,:) = 0.0
qtet2(:) = 0.0
frame: do i = 1, num_frames
        call cpu_time(t1)
        r1(:,:) = 0.0
        r1(:,:) = r(i,:,:)
        r2(:,:,:) = 0.0
        qtet_sum = 0.0

        outer: do j = 1, num_mol
                        inner: do k = j, num_mol
                                        !write(*,*) j,k
                                        ! we want distance between all the atoms for 1 frame
                                        r2(j,k,:) = r1(k,:) - r1(j,:)
                                        !write(300,*) j,k,r2(j,k,:) 
                                        ! MIC
                                        r2(j,k,:) = r2(j,k,:) - (box_len(i) * anint(r2(j,k,:)/box_len(i)))
                                       ! write(300,*) j,k,r2(j,k,:)
                                       ! distance calculation and saved
                                        r2(k,j,:) = -r2(j,k,:)
                                        dist(j,k) = sqrt(sum(r2(j,k,:)**2))
                                        dist(k,j) = dist(j,k)
                                     !   write(200,*) i,j,k,dist(j,k) 
                        end do inner
                        !do l = 1,num_mol
                        !        write(200,*) i,j,l,dist(j,l)
                        !end do
                        min_index = 0
                        temp = 0.0
                        temp_r2(:) = 0.0
                        do l = 1,5 ! to sort only first five distances
                                min_index = l
                                do k = l, num_mol
                                        if (dist(j,k) .lt. dist(j,min_index)) then
                                                min_index = k
                                        end if
                                end do
                                ! now swap to 1st value
                                temp = dist(j,min_index) !assign temp = minimum value
                                dist(j,min_index) = dist(j,l) ! assign the 8th place the value of 1st place
                                dist(j,l) = temp ! assign the temp value to 1st place
                                ! similarly swapping coordinates
                                temp_r2(:) = r2(j, min_index, :)
                                r2(j, min_index, :) = r2(j, l, :)
                                r2(j, l, :) = temp_r2(:)
                        end do
                        !do l = 1,num_mol
                        !        write(200,*) i,j,l,dist(j,l)
                        !end do
                        ! calculating qtet
                        sum1 = 0.0
                        sum2 = 0.0
                        do l = 2,4
                                do m = l+1,5
                                        theta = dot_product(r2(j,l,:),r2(j,m,:))/(dist(j, l) * dist(j, m))
                                       ! write(22,*) theta
                                        sum2 = (theta + 0.33d0)**2.0d0
                                        !write(23,*) sum2
                                        sum1 = sum1 + sum2
                                        !write(24,*) sum1
                                 end do
                        end do  
                        qtet1(i,j) = 1.0d0 - ((0.375d0)*sum1)
                        !write(51,*) i, qtet1(i,j)
                        qtet_sum = qtet_sum + qtet1(i,j)
        end do outer

        qtet2(i) = qtet_sum/num_mol
         !write(55,*) i, qtet2(i)
        call cpu_time(t2)
        if (mod(i,1000) == 0) then
                write(*,*) 'Time left in minutes is', ((t2-t1)*(num_frames-i)/60) , 'minutes'
        endif
end do frame


!do i = 1, num_mol
!        do j = 1, num_frames
!                write(50, *)qtet1(j, i)
!        end do
!end do
!close(50)

q_tot = sum(qtet2)/num_frames
write(*,*) 'Average Qtet over all frames', q_tot
write(100,*)  q_tot
write(*, *) 'Output file is output_qtet.txt'
print*, "Code executes successfully, Well DONE"
end program qtet
