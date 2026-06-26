from __future__ import annotations

from oracle_rovib import ORACLE_XYZ_ROTATIONAL_SCHEMA, ORACLE_XYZ_VIBRATIONAL_SCHEMA
from oracle_thermo import ORACLE_XYZ_THERMO_SCHEMA


def test_rovib_and_thermo_schema_constants_are_oracle_xyz_sections():
    assert ORACLE_XYZ_ROTATIONAL_SCHEMA == "oracle.xyz.rotational.v1"
    assert ORACLE_XYZ_VIBRATIONAL_SCHEMA == "oracle.xyz.vibrational.v1"
    assert ORACLE_XYZ_THERMO_SCHEMA == "oracle.xyz.thermo.v1"
