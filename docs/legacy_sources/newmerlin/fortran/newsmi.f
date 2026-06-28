      Implicit Integer (A-Z) 
      Parameter(MxAt=100,MxBnd=4,MxRng=10,MxBrn=10)
      Dimension IAn(MxAt),IArom(MxAt),IRing(MxAt),IBrnch(MxAt)
      Dimension NPi(MxAt),NBond(MxAt),IBond(MxBnd,MxAt)
      Dimension IniAt(MxAt),IniRng(MxRng),IEnRng(MxRng)
      Dimension IniBrn(MxBrn),IEnBrn(MxBrn)
      Character*100 SMILES
      Character*2 AtSymb(MxAt)
      Real*8 Charge(MxAt),BndOrd(MxBnd,MxAt)
C Define I/O files
      In=5
      Iout=6
C
      Read(In,*) SMILES
      call SMIF77(IOut,SMILES,MxBnd,NAtoms,NRing,NBrnch,AtSymb,IAn,
     $  IArom,NPi,IRing,IBrnch,NBond,IBond,IniAt,IniRng,IEnRng,
     $  IniBrn,IEnBrn,Charge,BndOrd)
      End
*Deck SMIF77
      SUBROUTINE SMIF77(Iout,SMILES,MxBnd,NAtoms,NRing,NBrnch,AtSymb,
     $  IAn,IArom,NPi,IRing,IBrnch,NBond,IBond,IniAt,IniRng,IEnRng,
     $  IniBrn,IEnBrn,Charge,BndOrd)
      IMPLICIT Integer (A-N) 
      Dimension IAn(*),NPi(*),IArom(*),IRing(*),IBrnch(*),NBond(*)
      Dimension IniAt(*),IBond(MxBnd,*)
      Dimension IniRng(*),IEnRng(*)
      Dimension IniBrn(*),IEnBrn(*)
      Character*(*) SMILES
      Character*2 AtSymb(*)
      Character*1 Star
      Real*8 Charge(*),BndOrd(MxBnd,*)
C
      IOut=6      
C Parse the SMILES and determine Atom properties
      CALL ParSMI(IOut,SMILES,NAtoms,NBrack,NRing,NBrnch,IArom,IniAt,
     $  Charge,AtSymb)
C Assign Brackets
C     CALL AssSqu(IOut,SMILES,NAtoms,IArom,Charge,AtSymb)
      write(IOut,'('' The Molecule contains'',I3,'' Atoms,'',I2,
     $ '' Rings'',I2,'' Branches and'',I2,'' Brackets'')') NAtoms,
     $ NRing,NBrnch,NBrack
C Assign Atomic Numbers
      Do IAt=1,NAtoms
       CALL FilIAn(AtSymb(IAt),IAn(IAt))
       write(IOut,'(A2,I3,I2)') AtSymb(IAt),IAt,IAn(IAt)
      end do
C Assign Rings
      CALL AssRng(IOut,SMILES,NAtoms,NRing,IniAt,IRing,IniRng,IEnRng)
C Assign Branches
      CALL AssBrn(IOut,SMILES,NAtoms,NBrnch,IniAt,IBrnch,IniBrn,IEnBrn) 
      write(IOut,'('' '')')
C Build connectivity
      CALL BldCon(IOut,SMILES,MxBnd,NAtoms,NPi,IArom,NRing,NBrnch,IniAt,
     $  IRing,IniRng,IEnRng,IBrnch,IniBrn,IEnBrn,NBond,IBond,
     $  BndOrd,AtSymb)
C Add Hydrogen Atoms
      NHeavy=NAtoms
      CALL AddHyd(IOut,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,IBond,BndOrd,
     $  AtSymb)
C Print Results
      CALL SmiRes(IOut,MxBnd,NAtoms,NPi,IArom,IBrnch,IRing,NBond,IBond,
     $  BndOrd,AtSymb)
      RETURN
      END
*Deck AssSqu
      Subroutine AssSqu(IOut,SMILES,NAtoms,IArom,Charge,AtSymb)
      Implicit Integer (A-Z)
      Dimension IArom(*)
      Real*8 Charge(*) 
      Character*(*) SMILES
      Character*1 El2
      Character*2 El,AtSymb(*)
      Logical InSqu
      LenSM=LEN_TRIM(SMILES)
      InSqu=.false.
      Do ISM=1,LenSM
       if(SMILES(ISM:ISM).eq.'[') then
        InSqu=.true.
        Cycle
       ElseIf(SMILES(ISM:ISM).eq.']') then
        InSqu=.false.
        Cycle
       End If
       If(InSqu) then
        If(SMILES(ISM:ISM).eq.'H') then
         nH=1
         El(1:1)='H' 
         El2=SMILES((ISM+1):(ISM+1))
         if(El2.ge.'1'.and.El2.le.'9') nH=IChar(El2)-IChar('0')
         do iH=1,nH
          NAtoms=NAtoms+iH
          IArom(NAtoms)=0
          Charge(NAtoms)=0.0D0
          AtSymb(NAtoms)='H '         
         end do
        End If
       End If
      End Do
      Return
      End 
*Deck ParSMI
      SUBROUTINE ParSMI(IOut,SMILES,NAtoms,NBrack,NRing,NBrnch,IArom,
     $  IniAt,Charge,AtSymb)
      Implicit Integer (A-Z)
      Dimension IArom(*),IniAt(*)
      Character*(*) SMILES
      Character*2 El,AtSymb(*)
      Real*8 Charge(*)
      Logical InSqu
      NAtoms=0
      NBrack=0
      NRing=0
      NBrnch=0
      IFndR=0
      IRing=0
      InSqu=.false.
      L=LEN_TRIM(SMILES)
      SMILES(L+1:L+1)=' '
      DO 10 I=1,L
       If(SMILES(I:I).eq.'[') then
        InSqu=.true.
        go to 10
       ElseIf(SMILES(I:I).eq.']') then
        InSqu=.false.
        goto 10
       End If
       If(InSqu) go to 10
       El(1:1)=SMILES(I:I)
       El(2:2)=SMILES(I+1:I+1)
       if(El(1:1).eq.'c') then 
        El(1:1)='C'
        If(El(2:2).ge.'1'.and.El(2:2).le.'9') then
         IRing=Ichar(El(2:2))-IChar('0')
         If(IRing.gt.IfndR) IFndR=IRing         
        EndIf
        El(2:2)=' '
        NAtoms=NAtoms+1
        IArom(NAtoms)=1
       elseif (El(1:1).eq.'n') then
        El(1:1)='N'
        If(El(2:2).ge.'1'.and.El(2:2).le.'9') then
         IRing=Ichar(El(2:2))-IChar('0')
         If(IRing.gt.IfndR) IFndR=IRing
        EndIf
        El(2:2)=' '
        NAtoms=NAtoms+1
        IArom(NAtoms)=1
       else
        if(El(1:1).ge.'a'.and.El(1:1).le.'z') goto 10
        if(El(1:1).ge.'1'.and.El(1:1).le.'9') goto 10 
        if(El(1:1).eq.'('.or.El(1:1).eq.')') goto 10 
        if(El(1:1).eq.'['.or.El(1:1).eq.']') goto 10  
        if(El(1:1).eq.'+'.or.El(1:1).eq.'-') goto 10  
        if(El(1:1).eq.'-'.or.El(1:1).eq.'='.or.El(1:1).eq.'#') goto 10  
        NAtoms=NAtoms+1
        Charge(NAtoms)=0.0D0
        IArom(NAtoms)=0 
        if(El(2:2).ge.'1'.and.El(2:2).le.'9') then
         IRing=Ichar(El(2:2))-IChar('0') 
         If(IRing.gt.IfndR) IFndR=IRing     
         El(2:2)=' '
        elseif(El(2:2).eq.'+') then 
         Charge(NAtoms)=1.0d0
         El(2:2)=' '
        elseif(El(2:2).eq.'-') then
         Charge(NAtoms)=-1.0d0
         El(2:2)=' '
        endif
        if(El(2:2).eq.'c'.or.El(2:2).eq.'n') El(2:2)=' '
        if(El(2:2).eq.'('.or.El(2:2).eq.')') El(2:2)=' '
        if(El(2:2).eq.'['.or.El(2:2).eq.']') El(2:2)=' '
        if(El(2:2).eq.'='.or.El(2:2).eq.'#') El(2:2)=' '
        if(El(2:2).ge.'A'.and.El(2:2).le.'Z') El(2:2)=' '
       endif 
       AtSymb(NAtoms)=El
       IniAt(NAtoms)=I
   10 continue   
      NRing=IfndR
      Do 20 I=1,L
       If(SMILES(I:I).eq.'(') NBrnch=NBrnch+1      
       If(SMILES(I:I).eq.'[') NBrack=NBrack+1
   20 continue
      return
      end
*Deck AssRng
      Subroutine AssRng(IOut,SMILES,NAtoms,NRing,IniAt,IRing,IniRng,
     $  IEnRng)
      Implicit Integer (A-Z)
      Character*(*) SMILES
      Logical Open
      Dimension IniAt(*),IRing(*),IniRng(*),IEnRng(*)
C Initialize all atoms to not be part of any ring
      Do IAt = 1, NAtoms
       IRing(IAt)=0
      end do
C Initialize ring  numbers
      Do i=1,NRing
       IniRng(i)=0
       IEnRng(i)=0
      end do
      NFRng=0
C Traverse the SMILES string to assign rings
      L=LEN_TRIM(SMILES)
      Open=.false.
      Do i=1,L
       if(SMILES(i:i).ge.'1'.and.SMILES(i:i).le.'9') then
        if(.not.Open) then
         NFRng=NFRng+1
         IniRng(NFRng)=i-1
         Open=.true.
        Else
         IEnRng(NFRng)=i-1
         Open=.false.
        end if
       end if
      end do
C Check the number of Rings
      If(NFRng.ne.NRing) then
       write(IOut,'(I2,'' Rings found instead of the expected'',I2)')
     $   NFRng,NRing
       stop
      endif
C Assign atoms to rings
      Do IAt=1,NAtoms
       Do j=1,NRing
        Ini=IniRng(j)
        IEnd=IEnRng(j)
        If(IniAt(IAt).ge.Ini.and.IniAt(IAt).le.IEnd) IRing(IAt)=j 
       End Do
      End Do
      Return
      End
*Deck AssBrn
      Subroutine AssBrn(IOut,SMILES,NAtoms,NBrnch,IniAt,IBrnch,IniBrn,
     $  IEnBrn)
      Implicit Integer (A-Z)
      Character*(*) SMILES
      Logical Open
      Dimension IniAt(*),IBrnch(*),IniBrn(*),IEnBrn(*)
C Initialize all atoms to not be part of any branch
      Do IAt = 1, NAtoms
       IBrnch(IAt)=0
      end do
C Initialize branch numbers
      Do i=1,NBrnch
       IniBrn(i)=0
       IEnBrn(i)=0
      end do
      NFBrn=0
C Traverse the SMILES string to assign branches
      LenSM=LEN_TRIM(SMILES)
      Open=.false.
      Do ISM=1,LenSM
       if(SMILES(ISM:ISM).eq.'(') then
         NFBrn=NFBrn+1
         IniBrn(NFBrn)=ISM+1
       ElseIf(SMILES(ISM:ISM).eq.')') then
        IEnBrn(NFBrn)=ISM-1
       end if
      end do
C Check the number of Branches
      If(NFBrn.ne.NBrnch) then
       write(IOut,'(I2,'' Branches found instead of the expected'',I2)')
     $   NFBrn,NBrnch
       stop
      endif
C Assign atoms to Branches
      Do IAt=1,NAtoms
       Do JBr=1,NBrnch
        Ini=IniBrn(JBr)
        IEnd=IEnBrn(JBr)
        If(IniAt(IAt).ge.Ini.and.IniAt(IAt).le.IEnd) IBrnch(IAt)=JBr
       End Do
      End Do
      Return
      End
*Deck BldCon
      Subroutine BldCon(IOut,SMILES,MxBnd,NAtoms,NPi,IArom,NRing,NBrnch,
     $  IniAt,IRing,IniRng,IEnRng,IBrnch,IniBrn,IEnBrn,NBond,IBond,
     $  BndOrd,AtSymb)
      Implicit Integer (A-Z)
      Dimension IniAt(*),IRing(*),IniRng(*),IEnRng(*)
      Dimension IBrnch(*),IniBrn(*),IEnBrn(*)
      Dimension NPi(*),IArom(*),NBond(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
      Character*(*) SMILES
      Character*2 AtSymb(*)
      Logical FndBnd
      NPi(NAtoms)=0
      do 10 IAt=1,NAtoms-1
       NPi(IAt)=0
       JAt=IAt+1       
       If(IBrnch(JAt).lt.IBrnch(IAt)) go to 10 
       InI=IniAt(IAt)
       IEnd=IniAt(JAt)
       NBond(IAt)=NBond(IAt)+1
       NBI=NBond(IAt)
       NBond(JAt)=NBond(JAt)+1
       NBJ=NBond(JAt)
       IBond(NBI,IAt)=JAt
       IBond(NBJ,JAt)=IAt
       BndOrd(NBI,IAt)=1.0d0
       BndOrd(NBJ,JAt)=1.0d0
       If(IArom(IAt).ne.0.and.IArom(JAt).ne.0) then
        BndOrd(NBI,IAt)=1.5d0
        BndOrd(NBJ,JAt)=1.5d0        
        NPi(IAt)=1
        NPi(JAt)=1
       EndIf
       If((IEnd-Ini).eq.2) then
        If(SMILES((IEnd-1):(IEnd-1)).eq.')') goto 10
        If(SMILES((IEnd-1):(IEnd-1)).eq.'=') then
         BndOrd(NBI,IAt)=2.0d0
         BndOrd(NBJ,JAt)=2.0d0
         NPi(IAt)=NPi(IAt)+1
         NPi(JAt)=NPi(IAt)+1
        ElseIf(SMILES(IEnd-1:IEnd-1).eq.'#') then
         BndOrd(NBI,IAt)=3.0d0
         BndOrd(NBJ,JAt)=3.0d0
         NPi(IAt)=2
         NPi(JAt)=2
        EndIf
       EndIf
       If((IEnd-Ini).eq.3) then
        If(SMILES((IEnd-2):(IEnd-2)).ne.'(') go to 10
        If(SMILES((IEnd-1):(IEnd-1)).eq.'=') then
         BndOrd(NBI,IAt)=2.0d0
         BndOrd(NBJ,JAt)=2.0d0
         NPi(IAt) = 1
         NPi(JAt) = 1          
        ElseIf(SMILES(IEnd-1:IEnd-1).eq.'#') then
         BndOrd(NBI,IAt)=3.0d0
         BndOrd(NBJ,JAt)=3.0d0
         NPi(IAt)=2
         NPi(JAt)=2
        EndIf
       EndIf
   10 continue
C Set bonds closing rings
      If(NRing.eq.0) go to 20
      Do IRefRn=1,NRing
       call ClRing(MxBnd,NAtoms,IRefRn,IArom,IRing,IniAt,IniRng,IEnRng,
     $  NBond,IBond,BndOrd) 
      End Do
C Set bonds before and after branches
   20 If(NBrnch.eq.0) return
      call BrnRef(MxBnd,NAtoms,NBrnch,NPi,IArom,IBrnch,NBond,IBond,
     $  BndOrd)
      Return
      End
*Deck BrnRef
      Subroutine BrnRef(MxBnd,NAtoms,NBrnch,NPi,IArom,IBrnch,NBond,
     $  IBond,BndOrd)
      Implicit Integer (A-Z)
      Dimension NPi(*),IArom(*),IBrnch(*),NBond(*),IBond(MxBnd,*) 
      Real*8 BndOrd(MxBnd,*)
      Logical Diff1,Fnd
      do 10 IRef=1,NBrnch
       Diff1=.false.
       Fnd=.false.
       do 20 IAt=1,NAtoms
        If(Fnd) go to 20
        If(IBrnch(IAt).eq.IRefi.and..not.Diff1) then
         Diff1=.true.
         IAt1=IAt-1
         go to 20
        ElseIf(Diff1.and.IBrnch(IAt).eq.(IRef-1)) then
         IAt2=IAt
         Fnd=.true.
        EndIf
   20  Continue
       NBond(IAt1)=NBond(IAt1)+1
       NBond(IAt2)=NBond(IAt2)+1
       IBond(NBond(IAt1),IAt1)=IAt2
       IBond(NBond(IAt2),IAt2)=IAt1
       BndOrd(NBond(IAt1),IAt1)=1.0d0 
       BndOrd(NBond(IAt2),IAt2)=1.0d0 
       If(IArom(IAt1).eq.1.and.IArom(IAt2).eq.1) then
        BndOrd(NBond(IAt1),IAt1)=1.5d0
        BndOrd(NBond(IAt2),IAt2)=1.5d0
       EndIf
   10 continue
      Return
      End
*Deck SmiRes
      Subroutine SmiRes(IOut,MxBnd,NAtoms,NPi,IArom,IBranch,IRing,
     $  NBond,IBond,BndOrd,AtSymb)
      Implicit Real*8 (A-H,O-Z)
      Dimension NPi(*),IArom(*),NBond(*),IBond(MxBnd,*)
      Dimension IBranch(*),IRing(*)
      Dimension BndOrd(MxBnd,*)
      Character*2 AtSymb(*)
      Character Star*1
      write(IOut,'(/,'' Aromatic Atoms are marked by * '')')
      write(IOut,'('' Atom Branch Ring  Pi El.'',4X,''Bonded Atoms'',
     $  10X,''Bond Orders'')')
      DO IAt = 1, NAtoms
       NBI = NBond(IAt)
       If(IArom(IAt).ne.0) then
        star='*'
       Else
        star=' '
       EndIf
       If(NBI.eq.1) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,2X,I2,7X,I3,12X,F8.4)')
     $     AtSymb(IAt),IAt,IBranch(IAt),IRing(IAt),Star,NPi(IAt),
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.2) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,2X,I2,7X,2I3,9X,2F8.4)')
     $     AtSymb(IAt),IAt,IBranch(IAt),IRing(IAt),Star,NPi(IAt),
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.3) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,2X,I2,7X,3I3,6X,3F8.4)')
     $     AtSymb(IAt),IAt,IBranch(IAt),IRing(IAt),Star,NPi(IAt),
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.4) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,2X,I2,7X,4I3,3X,4F8.4)')
     $     AtSymb(IAt),IAt,IBranch(IAt),IRing(IAt),Star,NPi(IAt),
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       Else
         write(IOut,'('' Wrong Number of Bonds ('',I1,'' ) for Atom'',
     $     I3)') NBI,IAt
          stop
         EndIf
      END DO
      write(IOut,'('' '')')
      Return
      End
*Deck AddHyd
      Subroutine AddHyd(IOut,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,
     $  IBond,BndOrd,AtSymb)
      Implicit Integer (A-Z) 
      Dimension NBond(*),IAn(*),NPi(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
      Character*2 AtSymb(*)
      IHAt = NHeavy
      DO 10 IAt=1,NHeavy
       MissH=NFree(IOut,IAt,IAn(IAt),NBond(IAt),NPi(IAt))  
       IF(MissH.gt.0) then 
        Do 20 JAt = 1, MissH
         IHAt = IHAt + 1
         IAn(IHAt)=1
         AtSymb(IHAt)='H '
         NBond(IAt)=NBond(IAt)+1
         NBond(IHAt)=1
         IBond(NBond(IAt),IAt)=IHAt 
         BndOrd(NBond(IAt),IAt)=1.0d0
         IBond(1,IHAt)=IAt
         BndOrd(1,IHAt)=1.0d0
  20    Continue
       EndIf 
  10  Continue
      NAtoms=IHAt
      Return
      End
*Deck NFree
      Integer Function NFree(IOut,IAt,IZAt,NSigma,NPi)
      Implicit Integer (A-N)
      IAv=IZAt
      If(IZAt.gt.18) then
       write(Iout,'('' Atom'',I3,'' has a not allowed Atomic Number:'',
     $   I2)')IAt,IZAt 
       Stop
      ElseIf(IZat.gt.10) then
       NAv=18-IZAt
      ElseIf(IZAt.gt.2) then
       NAv=10-IZAt
      EndIf
      NFree=NAv-NSigma-NPi
      If(NFree.lt.0) then
       write(IOut,'('' Atom'',I3,'' with Atomic Number'',I2,'' has'',I2,
     $ '' sigma and'',I2,'' pi bonds'')') IAt,IZAt,NSigma,NPi  
       Stop
      EndIf
      Return
      End
C ---------------------------------------------------------------
C                   Routines from Tools1
C ---------------------------------------------------------------
*Deck FilIAn
      Subroutine FilIAn(El,IAn)
      Integer El2IAn
      Character*10 Num
      Character*2 El 
      Save Num
      Data Num/'0123456789'/
  100 Format(A1)
      IAn=El2IAN(.True.,El)
      If(IAn.gt.-1) Return
      If(El(2:2).eq.' ') then
       do 10 i=1,10
        If(El(1:1).eq.Num(i:i)) IAn=IAn+i
   10  continue
      Else
       Do 20 i=1,10
        If(El(2:2).eq.Num(i:i)) IAn=Ian+i
        If(El(1:1).eq.Num(i:i)) IAn=IAn+(i-1)*10
   20  Continue 
      EndIf
      Return
      End
*Deck El2IAN
      Integer Function El2IAN(NoCase,AtSymb)
      Implicit Integer(A-Z)
C
C     Returns the atomic number for a given atomic symbol.
C
C Input:
C     NoCase : Ignore case while looking for atomic number
C     AtSymb : Atomic symbol (character)
C     
C     Input
      Character AtSymb*2
      Logical NoCase
C     Local
      Integer i
      Character AllSmb*222, ASymb*2
C     
      If(NoCase) then
        AllSmb = ' ' // 'H HELIBEB C N O F NENAMGALSIP S CLARK CASCTIV '
     $  //'CRMNFECONICUZNGAGEASSEBRKRRBSRY ZRNBMOTCRURHPDAGCDINSNSBTEI '
     $  //'XECSBALACEPRNDPMSMEUGDTBDYHOERTMYBLUHFTAW REOSIRPTAUHGTLPBBI'
     $  //'POATRNFRRAACTHPAU NPPUAMCMBKCFESFMMDNOLRRFDBSGBHHSMTDS'
        Call LinUpC(AtSymb,ASymb)
      else
        AllSmb = ' ' // 'H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV '
     $  //'CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI '
     $  //'XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBi'
     $  //'PoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDs'
        ASymb = AtSymb
      endIf
      i = Index(AllSmb,ASymb)
      If(i.eq.0) then
        El2IAN = -1
      else
        El2IAN = i/2
      endIf
      Return
      End
*Deck LinUpC 
      Subroutine LinUpc(InStr,OutStr)
      Implicit Integer(A-Z)
C    
C     Translate a character string to upper case.
C    
      Character*(*) InStr, OutStr
C     
      IUA = IChar('A')
      ILA = IChar('a')
      ILZ = IChar('z')
      LenO = Min(Len(InStr),Len(OutStr))
      Do 10 I = 1, LenO
        ICI = IChar(InStr(I:I))
        If(ICI.ge.ILA.and.ICI.le.ILZ) ICI = ICI + IUA - ILA
        OutStr(I:I) = Char(ICI)
   10 Continue
      If(LenO.lt.Len(OutStr)) OutStr(LenO+1:) = ' '
      Return
      End
*Deck ClRing
      Subroutine ClRing(MxBnd,NAtoms,IRefRn,IArom,IRing,IniAt,IniRng,
     $  IEnRng,NBond,IBond,BndOrd)
      Implicit Integer (A-Z)
      Dimension IArom(*),IRing(*),IniAt(*),IniRng(*),IEnRng(*)
      Dimension NBond(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
      Do 10 IAt=1,NAtoms-1
       IAtIni=IniAt(IAt)
       IRI=IRing(IAt)
       If(IRI.ne.IRefRn) goto 10
       INRI=IniRng(IRI)
       IERI=IEnRng(IRI)
       If(IAtIni.ne.INRI) goto 10
       Do 20 JAt=IAt+1,NAtoms
        JAtIni=IniAt(JAt)
        If(JAtIni.ne.IERI) goto 20
        NBond(IAt)=NBond(IAt)+1
        NBond(JAt)=NBond(JAt)+1
        IBond(NBond(IAt),IAt)=JAt
        IBond(NBond(JAt),JAt)=IAt
        If(IArom(IAt).ne.0.and.IArom(JAt).ne.0) then
         BndOrd(NBond(IAt),IAt)=1.5D0
         BndOrd(NBond(JAt),JAt)=1.5D0
        Else
         BndOrd(NBond(IAt),IAt)=1.0D0
         BndOrd(NBond(JAt),JAt)=1.0D0
        EndIf
   20  continue
   10 continue
      Return
      End

