from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


def test_fortran_fragment_tric_bmat_compiles_and_runs(tmp_path):
    gfortran = shutil.which("gfortran")
    if gfortran is None:
        pytest.skip("gfortran is not available")

    root = Path(__file__).resolve().parents[1]
    source = root / "engines" / "fortran" / "gicforge" / "frag_tric_bmat.f"
    driver = tmp_path / "driver.f"
    executable = tmp_path / "driver"
    driver.write_text(
        """
      Program TFRAG
      Implicit Real*8(A-H,O-Z)
      Integer FAT(3),RAT(3)
      Dimension C(3,6),B(18)
      Data FAT /4,5,6/
      Data RAT /1,2,3/
C
      Do 10 I=1,18
       B(I)=0.0D0
   10 Continue
      Do 20 I=1,6
       Do 30 J=1,3
        C(J,I)=0.0D0
   30  Continue
   20 Continue
C
      C(1,1)=0.00D0
      C(2,1)=0.00D0
      C(3,1)=0.00D0
      C(1,2)=0.96D0
      C(2,2)=0.00D0
      C(3,2)=0.00D0
      C(1,3)=-0.24D0
      C(2,3)=0.93D0
      C(3,3)=0.00D0
      C(1,4)=0.00D0
      C(2,4)=0.00D0
      C(3,4)=3.20D0
      C(1,5)=0.96D0
      C(2,5)=0.00D0
      C(3,5)=3.20D0
      C(1,6)=-0.24D0
      C(2,6)=0.93D0
      C(3,6)=3.20D0
C
      Call ORCFTRN(6,3,FAT,3,RAT,1,B)
      If(DAbs(B(10)-1.0D0/3.0D0).gt.1.0D-10) Stop 1
      If(DAbs(B(1)+1.0D0/3.0D0).gt.1.0D-10) Stop 2
C
      Call ORCFCDI(6,3,FAT,3,RAT,C,B,VAL,IFAIL)
      If(IFAIL.ne.0) Stop 3
      If(DAbs(VAL-3.20D0).gt.1.0D-10) Stop 4
      If(DAbs(B(12)-1.0D0/3.0D0).gt.1.0D-10) Stop 5
      If(DAbs(B(3)+1.0D0/3.0D0).gt.1.0D-10) Stop 6
C
      Call ORCFROT(6,3,FAT,3,RAT,5,6,2,3,1,C,B,VAL,IFAIL)
      If(IFAIL.ne.0) Stop 7
      If(DAbs(VAL).gt.1.0D-10) Stop 8
      SUM=0.0D0
      Do 40 I=1,18
       If(B(I).ne.B(I)) Stop 9
       SUM=SUM+DAbs(B(I))
   40 Continue
      If(SUM.le.1.0D-8) Stop 10
      End
""",
        encoding="ascii",
    )

    subprocess.run(
        [
            gfortran,
            "-std=legacy",
            "-Wall",
            "-Wextra",
            "-fcheck=all",
            str(source),
            str(driver),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=root,
    )
    subprocess.run([str(executable)], check=True, cwd=root)
