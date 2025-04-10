module mo_math
    implicit none

    interface swap
        module procedure :: swap_real8, swap_integer
    end interface

    interface atan2
        module procedure :: atan2_complex
    end interface

    interface qsort
        module procedure :: qsort_real8, qsort_integer
    end interface

contains

    subroutine init_rand(seed)
        implicit none

        ! para list
        integer :: seed

        ! local
        integer :: n
        integer, allocatable, dimension(:) :: seed_array

        call random_seed(size=n)
        allocate( seed_array(n) )
        seed_array = seed
        call random_seed(put=seed_array)
    end subroutine

    function randperm(n) result(redata)
        implicit none

        ! para list
        integer, intent(in) :: n

        ! result
        integer, dimension(n) :: redata

        ! local
        integer :: i, itemp
        real(8) :: rtemp

        redata(1) = 1
        do i=2, n
            call random_number(rtemp)
            itemp = floor(rtemp*i) + 1
            if ( i /= itemp ) then
                redata(i) = redata(itemp)
                redata(itemp) = i
            else
                redata(i) = i
            end if
        end do
    end function

    function randuvec(n) result(uvec)
        implicit none

        ! para list
        integer, intent(in) :: n

        ! result
        real(8) :: uvec(n)

        call random_number(uvec)
        uvec = 2 * uvec - 1.d0
        do while ( norm2(uvec) > 1.d0 )
            call random_number(uvec)
            uvec = 2 * uvec - 1.d0
        end do

        uvec = uvec / norm2(uvec)
    end function

    ! swap >
        subroutine swap_real8(x, y)
            implicit none

            ! para list
            real(8), intent(inout) :: x, y

            ! local
            real(8) :: tmp

            tmp = x
            x   = y
            y   = tmp
        end subroutine

        subroutine swap_integer(x, y)
            implicit none

            ! para list
            integer, intent(inout) :: x, y

            ! local
            integer :: tmp

            tmp = x
            x   = y
            y   = tmp
        end subroutine
    ! swap <

    pure function mean(a) result(re)
        implicit none

        ! para list
        real(8), intent(in), dimension(:) :: a

        ! result
        real(8) :: re

        re = sum(a) / size(a)
    end function

    pure function std(a) result(re)
        implicit none

        ! para list
        real(8), intent(in), dimension(:) :: a

        ! result
        real(8) :: re

        integer :: ilen
        real(8) :: mean_of_a

        ilen = size(a)

        mean_of_a = sum(a) / ilen

        re = sum( (a-mean_of_a)**2 ) / ilen
        re = sqrt(re)
    end function

    pure function corr(a, b) result(re)
        implicit none

        ! para list
        real(8), intent(in), dimension(:) :: a, b

        ! result
        real(8) :: re

        ! local
        integer :: ilen
        real(8) :: sigma_of_a, sigma_of_b

        ilen = size(a)

        sigma_of_a = std(a)
        sigma_of_b = std(b)

        re = sum(a*b)/ilen - mean(a)*mean(b)
        re = re / std(a) / std(b)
    end function

    pure function dot(a,b) result(re)
        implicit none

        ! para list
        real(8), intent(in), dimension(:) :: a, b

        ! result
        real(8) :: re

        re = sum(a*b)
    end function

    pure function times2(a,b) result(re)
        implicit none

        ! para list
        real(8), intent(in)  :: a(2), b(2)

        ! result
        real(8) :: re

        re = a(1) * b(2) - a(2) * b(1)
    end function

    pure function times3(a,b) result(re)
        implicit none

        ! para list
        real(8), intent(in)  :: a(3), b(3)

        ! result
        real(8) :: re(3)

        re(1) = a(2) * b(3) - a(3) * b(2)
        re(2) = a(3) * b(1) - a(1) * b(3)
        re(3) = a(1) * b(2) - a(2) * b(1)
    end function

    pure function unitv(vector) result(uvector)
        implicit none

        ! para list
        real(8), intent(in),  dimension(:) :: vector

        ! result
        real(8), allocatable, dimension(:) :: uvector

        uvector = vector / norm2(vector)
    end function

    function sortperm(n, data) result(index)
        implicit none
        !===================================================================
        !
        !     sortrx -- sort, real input, index output
        !
        !
        !     input:  n     integer
        !             data  real
        !
        !     output: index integer (dimension n)
        !
        ! this routine performs an in-memory sort of the first n elements of
        ! array data, returning into array index the indices of elements of
        ! data arranged in ascending order.  thus,
        !
        !    data(index(1)) will be the smallest number in array data;
        !    data(index(n)) will be the largest number in data.
        !
        ! the original data is not physically rearranged.  the original order
        ! of equal input values is not necessarily preserved.
        !
        !===================================================================
        !
        ! sortrx uses a hybrid quicksort algorithm, based on several
        ! suggestions in knuth, volume 3, section 5.2.2.  in particular, the
        ! "pivot key" [my term] for dividing each subsequence is chosen to be
        ! the median of the first, last, and middle values of the subsequence;
        ! and the quicksort is cut off when a subsequence has 9 or fewer
        ! elements, and a straight insertion sort of the entire array is done
        ! at the end.  the result is comparable to a pure insertion sort for
        ! very short arrays, and very fast for very large arrays (of order 12
        ! micro-sec/element on the 3081k for arrays of 10k elements).  it is
        ! also not subject to the poor performance of the pure quicksort on
        ! partially ordered data.
        !
        ! created:  15 jul 1986  len moss
        !
        !===================================================================

        integer, intent(in) :: n
        real(8), intent(in) :: data(n)
        integer :: index(n)

        integer,dimension(31) :: lstk,rstk
        integer               :: istk
        integer               :: l,r,i,j,p,indexp,indext
        real(8)               :: datap

        ! quicksort cutoff
        !
        ! quit quicksort-ing when a subsequence contains m or fewer
        ! elements and finish off at end with straight insertion sort.
        ! according to knuth, v.3, the optimum value of m is around 9.

        integer, parameter  :: m = 9

        !===================================================================
        !
        ! make initial guess for index

        do i=1,n
            index(i)=i
        end do
        ! if array is short, skip quicksort and go directly to
        ! the straight insertion sort.

        if (n .le. m) goto 900

        !===================================================================
        !
        ! quicksort
        !
        ! the "qn:"s correspond roughly to steps in algorithm q,
        ! knuth, v.3, pp.116-117, modified to select the median
        ! of the first, last, and middle elements as the "pivot
        ! key" (in knuth's notation, "k").  also modified to leave
        ! data in place and produce an index array.  to simplify
        ! comments, let data[i]=data(index(i)).

        ! q1: initialize
        istk=0
        l=1
        r=n

        200 continue

        ! q2: sort the subsequence data[l]..data[r].
        !
        !     at this point, data[l] <= data[m] <= data[r] for all l < l,
        !     r > r, and l <= m <= r.  (first time through, there is no
        !     data for l < l or r > r.)

        i=l
        j=r

        ! q2.5: select pivot key
        !
        ! let the pivot, p, be the midpoint of this subsequence,
        ! p=(l+r)/2; then rearrange index(l), index(p), and index(r)
        ! so the corresponding data values are in increasing order.
        ! the pivot key, datap, is then data[p].

        p=(l+r)/2
        indexp=index(p)
        datap=data(indexp)

        if (data(index(l)) .gt. datap) then
            index(p)=index(l)
            index(l)=indexp
            indexp=index(p)
            datap=data(indexp)
        end if

        if (datap .gt. data(index(r))) then
            if (data(index(l)) .gt. data(index(r))) then
                index(p)=index(l)
                index(l)=index(r)
            else
                index(p)=index(r)
            end if
            index(r)=indexp
            indexp=index(p)
            datap=data(indexp)
        end if

        ! now we swap values between the right and left sides and/or
        ! move datap until all smaller values are on the left and all
        ! larger values are on the right.  neither the left or right
        ! side will be internally ordered yet; however, datap will be
        ! in its final position.

        300 continue

        ! q3: search for datum on left >= datap
        !
        ! at this point, data[l] <= datap.  we can therefore start scanning
        ! up from l, looking for a value >= datap (this scan is guaranteed
        ! to terminate since we initially placed datap near the middle of
        ! the subsequence).

        i=i+1
        if (data(index(i)).lt.datap) goto 300

        400 continue

        ! q4: search for datum on right <= datap
        !
        ! at this point, data[r] >= datap.  we can therefore start scanning
        ! down from r, looking for a value <= datap (this scan is guaranteed
        ! to terminate since we initially placed datap near the middle of
        ! the subsequence).

        j=j-1
        if (data(index(j)).gt.datap) goto 400

        ! q5: have the two scans collided?

        if (i.lt.j) then

        ! q6: no, interchange data[i] <--> data[j] and continue

            indext=index(i)
            index(i)=index(j)
            index(j)=indext
            goto 300
        else

        ! q7: yes, select next subsequence to sort
        !
        ! at this point, i >= j and data[l] <= data[i] == datap <= data[r],
        ! for all l <= l < i and j < r <= r.  if both subsequences are
        ! more than m elements long, push the longer one on the stack and
        ! go back to quicksort the shorter; if only one is more than m
        ! elements long, go back and quicksort it; otherwise, pop a
        ! subsequence off the stack and quicksort it.

        if (r-j .ge. i-l .and. i-l .gt. m) then
            istk=istk+1
            lstk(istk)=j+1
            rstk(istk)=r
            r=i-1
        else if (i-l .gt. r-j .and. r-j .gt. m) then
            istk=istk+1
            lstk(istk)=l
            rstk(istk)=i-1
            l=j+1
        else if (r-j .gt. m) then
            l=j+1
        else if (i-l .gt. m) then
            r=i-1
        else
        ! q8: pop the stack, or terminate quicksort if empty
            if (istk.lt.1) goto 900
                l=lstk(istk)
                r=rstk(istk)
                istk=istk-1
            end if
            goto 200
        end if

        900 continue

        !===================================================================
        !
        ! q9: straight insertion sort

        do 950 i=2,n
            if (data(index(i-1)) .gt. data(index(i))) then
                indexp=index(i)
                datap=data(indexp)
                p=i-1
        920     continue
                    index(p+1) = index(p)
                    p=p-1
                    if (p.gt.0) then
                        if (data(index(p)).gt.datap) goto 920
                    end if
                index(p+1) = indexp
            end if
        950 continue

        !===================================================================
        !
        !     all done
    end function

    ! qsort >
        ! for real
        recursive subroutine qsort_real8(data)
            implicit none

            ! para list
            real(8), intent(inout) :: data(:)

            ! local
            integer :: iq

            if ( size(data) > 1 ) then
                call partition_data(data, iq)
                call qsort(data(:iq-1))
                call qsort(data(iq:))
            end if
            contains
            subroutine partition_data(pdata, marker)
                implicit none

                ! para list
                real(8), intent(inout) :: pdata(:)
                integer, intent(out)   :: marker

                ! local
                integer :: i, j
                real(8) :: pivot

                i=0
                j=size(pdata)+1
                pivot = pdata((i+j)/2)

                do while (.true.)
                    i = i + 1
                    j = j - 1
                    do while ( pdata(i) < pivot )
                       i = i + 1
                    end do
                    do while ( pdata(j) >  pivot )
                        j = j - 1
                    end do
                    if ( i < j ) then
                        call swap(pdata(i),pdata(j))
                    else if ( i == j ) then
                        marker = i + 1; exit
                    else ! i>j
                        marker = i; exit
                    end if
                end do
            end subroutine
        end subroutine

        ! for integer
        recursive subroutine qsort_integer(data)
            implicit none

            ! para list
            integer, intent(inout) :: data(:)

            ! local
            integer :: iq

            if ( size(data) > 1 ) then
                call partition_data(data, iq)
                call qsort(data(:iq-1))
                call qsort(data(iq:))
            end if
            contains
            subroutine partition_data(pdata, marker)
                implicit none

                ! para list
                integer, intent(inout) :: pdata(:)
                integer, intent(out)   :: marker

                ! local
                integer :: i, j
                integer :: pivot

                i=0
                j=size(pdata)+1
                pivot = pdata((i+j)/2)

                do while (.true.)
                    i = i + 1
                    j = j - 1
                    do while ( pdata(i) < pivot )
                       i = i + 1
                    end do
                    do while ( pdata(j) >  pivot )
                        j = j - 1
                    end do
                    if ( i < j ) then
                        call swap(pdata(i),pdata(j))
                    else if ( i == j ) then
                        marker = i + 1; exit
                    else ! i>j
                        marker = i; exit
                    end if
                end do
            end subroutine
        end subroutine
    ! swap <

    pure function atan2_complex(x) result(re)
        implicit none

        ! para list
        complex(16), intent(in) :: x

        ! result
        real(8) :: re

        re = atan2(aimag(x),real(x))
    end function

    pure function vec_angle( vec1, vec2 ) result(angle)
        ! calc angle of vec1 and vec2
        implicit none

        ! para list
        real(8), dimension(2), intent(in) :: vec1, vec2

        ! result
        real(8) :: angle

        ! local
        complex(16) :: v1, v2

        v1 = cmplx(vec1(1),vec1(2))
        v2 = cmplx(vec2(1),vec2(2))

        angle = atan2( v2/v1 )
        !angle = atan2(aimag(v3),real(v3))
    end function

    pure function dvec_angle( vec1, vec2 ) result(angle)
        ! calc angle of vec1 and (vec2-vec1)
        implicit none

        ! para list
        real(8), dimension(2), intent(in) :: vec1, vec2

        ! result
        real(8) :: angle

        ! local
        real(8), dimension(2) :: vec3
        complex(16) :: v1, v2

        vec3 = vec2 - vec1

        v1 = cmplx(vec1(1),vec1(2))
        v2 = cmplx(vec3(1),vec3(2))

        angle = atan2( v2/v1 )

        !v3 = v2 / v1
        !angle = atan2(aimag(v3),real(v3))
    end function

    pure function findfirst(data,x) result(iloc)
        ! find first value in data equal to x
        ! return index
        implicit none

        ! para list
        integer, intent(in) :: data(:)
        integer, intent(in) :: x

        ! result
        integer :: iloc

        ! local
        integer :: n, istart, imed, iend, i

        n = size(data)
        iloc = 0

        istart = 1
        iend   = n
        imed   = (istart+n)/2

        do while ( iend - istart > 10 )
            if ( x <= data(imed) ) then
                iend = imed
                imed = (istart+iend)/2
            else
                istart = imed
                imed = (istart+iend)/2
            end if
        end do

        do i=istart, iend
            if ( data(i)==x ) then
                iloc = i
                exit
            end if
        end do
    end function

end module
