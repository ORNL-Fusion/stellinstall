      SUBROUTINE allocate_angles
      USE parambs
      USE vmec0
      IMPLICIT NONE
      INTEGER :: istat
      ALLOCATE(bfield(nthetah,nzetah), gsqrt_b(nthetah,nzetah),
     1  sinmi(-mbuse:mbuse,nthetah), sinnj(0:nbuse,nzetah),
     2  cosmi(-mbuse:mbuse,nthetah), cosnj(0:nbuse,nzetah),
     3  theta(nthetah), zetah(nzetah),
     4  sinmim(-mbuse:mbuse,nthetahm), sinnjm(0:nbuse,nzetahm),
     5  cosmim(-mbuse:mbuse,nthetahm), cosnjm(0:nbuse,nzetahm),
     6  b2obm(nthetah,nzetah), bfieldm(nthetahm,nzetahm), stat=istat)

      IF (istat .ne. 0) STOP 'allocation error in allocate_angles'

      END SUBROUTINE allocate_angles
