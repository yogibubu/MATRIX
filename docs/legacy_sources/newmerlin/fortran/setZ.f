      Dimension NBond(10),IZ(4,10),IBond(4,10)
      do 10 IAt=1,10
       NBond(IAt)=0
       do 20 ILeg=1,4
        IZ(ILeg,IAt)=0
        IBond(ILeg,IAt)=0
   20  continue    
   10 continue 
C Normale
      IZ(1,1)=0
      IZ(2,1)=0
      IZ(3,1)=0 
      IZ(1,2)=1
      IZ(2,2)=0
      IZ(3,2)=0
      IZ(1,3)=2
      IZ(2,3)=1
      IZ(3,3)=0
      NBond(1)=1
      NBond(2)=2
      NBond(3)=2
      IBond(1,1)=2
      IBond(1,2)=1
      IBond(2,2)=3
      IBond(1,3)=1
      IAt=4
      JAt=3       
      Call SetIZ(IAt,JAt,KAt,LAt,NBond,IBond,IZ)
      write(*,*) JAt,KAt,LAt,IZ(4,4)
      write(*,*) IZ(1,4),IZ(2,4),IZ(3,4),IZ(4,4)
C Branched
      IZ(1,1)=0
      IZ(2,1)=0
      IZ(3,1)=0
      IZ(1,2)=1
      IZ(2,2)=0
      IZ(3,2)=0
      IZ(1,3)=1
      IZ(2,3)=2
      IZ(3,3)=0
      NBond(1)=2
      NBond(2)=1
      NBond(3)=1
      IBond(1,1)=2
      IBond(2,1)=3
      IBond(1,2)=1
      IBond(1,3)=1
      IAt=4
      JAt=3
      Call SetIZ(IAt,JAt,KAt,LAt,NBond,IBond,IZ)
      write(*,*) JAt,KAt,LAt,IZ(4,4)
      write(*,*) IZ(1,4),IZ(2,4),IZ(3,4),IZ(4,4)
C Star
      IZ(1,1)=0
      IZ(2,1)=0
      IZ(3,1)=0
      IZ(1,2)=1
      IZ(2,2)=3
      IZ(3,2)=0
      IZ(1,3)=2
      IZ(2,3)=0
      IZ(3,3)=0
      NBond(1)=1
      NBond(2)=3
      NBond(3)=1
      IBond(1,1)=2
      IBond(2,1)=1
      IBond(2,2)=3
      IBond(1,3)=2
      IAt=4
      JAt=2
      Call SetIZ(IAt,JAt,KAt,LAt,NBond,IBond,IZ)
      write(*,*) JAt,KAt,LAt,IZ(4,4) 
      write(*,*) IZ(1,4),IZ(2,4),IZ(3,4),IZ(4,4)
      End
*Deck SetIZ
      Subroutine SetIZ(IAt,JAt,KAt,LAt,NBond,IBond,IZ)
      Dimension NBond(10),IZ(4,10),IBond(4,10)
      IZ(1,IAt)=JAt
      If(IAt.eq.2) return
      KAt=IZ(1,JAt)
      If(KAt.eq.0) then
       KAt=IBond(2,JAt)
       If(IBond(2,JAt).eq.IAt) KAt=IBond(1,JAt)
      EndIf
      IZ(2,IAt)=KAt 
      If(IAt.eq.3) return
C fourth atom
      If(KAt.eq.1) then
       If(NBond(KAt).gt.1) then
        LAt=IBond(1,KAt)
        If(LAt.eq.JAt) LAt=IBond(2,KAt)
        IZ(4,IAt)=0
       Else
        LAt=IBond(2,JAt)
        If(LAt.eq.KAt) LAt=IBond(1,JAt)
        IZ(4,IAt)=1
       EndIf
      Else
       If(NBond(KAt).gt.1) then
        LAt=IZ(1,KAt)
        IZ(4,IAt)=0
       Else
        LAt=IBond(2,JAt)
        If(LAt.eq.KAt) LAt=IBond(1,JAt)
        IZ(4,IAt)=1
       EndIf
      EndIf
      IZ(3,IAt)=LAt
      Return
      End
      
