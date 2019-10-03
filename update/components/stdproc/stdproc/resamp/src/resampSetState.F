!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! copyright: 2010 to the present, california institute of technology.
! all rights reserved. united states government sponsorship acknowledged.
! any commercial use must be negotiated with the office of technology transfer
! at the california institute of technology.
! 
! this software may be subject to u.s. export control laws. by accepting this
! software, the user agrees to comply with all applicable u.s. export laws and
! regulations. user has the responsibility to obtain export licenses,  or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.
! 
! installation and use of this software is restricted by a license agreement
! between the licensee and the california institute of technology. it is the
! user's responsibility to abide by the terms of the license agreement.
!
! Author: Giangi Sacco
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




        subroutine setStdWriter(varInt)
            use resampState
            implicit none
            integer*8 varInt
            ptStdWriter = varInt
        end

        subroutine setNumberFitCoefficients(varInt)
            use resampState
            implicit none
            integer varInt
            i_ma = varInt
        end

        subroutine setNumberRangeBin1(varInt)
            use resampState
            implicit none
            integer varInt
            npl = varInt
        end

        subroutine setNumberRangeBin2(varInt)
            use resampState
            implicit none
            integer varInt
            npl2 = varInt
        end

        subroutine setStartLine(varInt)
            use resampState
            implicit none
            integer varInt
            ist = varInt
        end

        subroutine setNumberLines(varInt)
            use resampState
            implicit none
            integer varInt
            nl = varInt
        end

        subroutine setNumberLinesImage2(varInt)
            use resampState
            implicit none
            integer varInt
            ndown2 = varInt
        end

        subroutine setFirstLineOffset(varInt)
            use resampState
            implicit none
            integer varInt
            istoff = varInt
        end

        subroutine setNumberRangeLooks(varInt)
            use resampState
            implicit none
            integer varInt
            NR = varInt
        end

        subroutine setNumberAzimuthLooks(varInt)
            use resampState
            implicit none
            integer varInt
            NAZ = varInt
        end

        subroutine setRadarWavelength(varInt)
            use resampState
            implicit none
            real*4 varInt
            WVL = varInt
        end

        subroutine setSlantRangePixelSpacing(varInt)
            use resampState
            implicit none
            real*4 varInt
            SLR = varInt
        end

        subroutine setFlattenWithOffsetFitFlag(varInt)
            use resampState
            implicit none
            integer varInt
            iflatten = varInt
        end

        subroutine setDopplerCentroidCoefficients(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                dopplerCoefficients(i) = array1d(i)
            enddo
        end

        subroutine setLocationAcross1(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_ranposV(i) = array1d(i)
            enddo
        end

        subroutine setLocationAcrossOffset1(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_ranoffV(i) = array1d(i)
            enddo
        end

        subroutine setLocationDown1(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_azposV(i) = array1d(i)
            enddo
        end

        subroutine setLocationDownOffset1(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_azoffV(i) = array1d(i)
            enddo
        end

        subroutine setSNR1(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_sigV(i) = array1d(i)
            enddo
        end

        subroutine setLocationAcross2(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_ranpos2V(i) = array1d(i)
            enddo
        end

        subroutine setLocationAcrossOffset2(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_ranoff2V(i) = array1d(i)
            enddo
        end

        subroutine setLocationDown2(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_azpos2V(i) = array1d(i)
            enddo
        end

        subroutine setLocationDownOffset2(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_azoff2V(i) = array1d(i)
            enddo
        end

        subroutine setSNR2(array1d,dim1)
            use resampState
            implicit none
            integer dim1,i
            double precision, dimension(dim1):: array1d
            do i = 1, dim1
                r_sig2V(i) = array1d(i)
            enddo
        end

