# Plate package for UEDGE, plates.py
# Extended from plate_d3d_10.py
# Changelog
# 200214 - file created


def plate_d3d_10():
    from uedge import grd
    # Import pointer to grd
 
    # for inboard half of mesh: (same as outboard half)
    grd.nplate1=5
    grd.gchange("grd.Mmod",0)

    grd.rplate1=[   1.600E+00,      1.27300E+00,    1.15310E+00,
                    1.01600E+00,    1.01600E+00
                ]

    grd.zplate1=[   2.34100E-01,    2.34100E-01,    2.34100E-01,  
                    3.71200E-01,    1.0000
                ]

    # for outboard half of mesh: (modified baffle/plate)
    grd.nplate2=10
    grd.gchange("grd.Mmod",0)

    grd.rplate2=[   2.13690E+00,    1.78570E+00,    1.76800E+00,  
                    1.76800E+00,    1.68100E+00,    1.67500E+00,  
                    1.67200E+00,    1.67200E+00,    1.55500E+00,  
                    1.21200E+00
                ]

    grd.zplate2=[   6.28600E-01,    4.25600E-01,    3.89300E-01,  
                    0.35,           0.35,           0.35, 
                    0.35,           0.35,           0.35, 
                    0.35
                ]

