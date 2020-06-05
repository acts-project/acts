
// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "TGeoBBox.h"
#include "TGeoBoolNode.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoNode.h"
#include "TGeoPgon.h"
#include "TGeoVolume.h"

namespace Acts {

namespace Test {

struct TGeoPixelBarrel {
  /// @brief this is a ROOT description of a silicon detector layer
  /// which is only used for unit testing
  ///
  /// @param volumeName The name of the search volume (top volume here)
  /// @param sensorName The name of the sensor to select
  ///
  /// The TopVolume is returned
  inline static TGeoVolume *construct(
      const std::string &volumeName = "PixelBarrel",
      const std::string &sensorName = "SiSensor") {
    new TGeoManager("TestTGeoGeometry", "TestTGeoGeometry");

    Double_t dx, dy, dz, phi1, dphi;
    Double_t par[20];
    Double_t rmin, rmax;
    Double_t a;
    Double_t thx, phx, thy, phy, thz, phz;
    Double_t z, density, radl, absl, w;
    Int_t nel, numed, nz, nedges;

    // MATERIALS, MIXTURES AND TRACKING MEDIA
    // Mixture: T_Air
    nel = 4;
    density = 0.001214;
    auto matAir = new TGeoMixture("T_Air", nel, density);
    a = 14.007;
    z = 7;
    w = 0.74939999999999996;  // N
    matAir->DefineElement(0, a, z, w);
    a = 15.999000000000001;
    z = 8;
    w = 0.2369;  // O
    matAir->DefineElement(1, a, z, w);
    a = 39.948;
    z = 18;
    w = 0.0129;  // AR
    matAir->DefineElement(2, a, z, w);
    a = 1.0079400000000001;
    z = 1;
    w = 0.0007999999;  // H
    matAir->DefineElement(3, a, z, w);
    matAir->SetIndex(30);
    // Medium: T_Air
    numed = 365;                   // medium number
    par[0] = 0;                    // isvol
    par[1] = -1;                   // ifield
    par[2] = 40;                   // fieldm
    par[3] = 20;                   // tmaxfd
    par[4] = 10000000000;          // stemax
    par[5] = 0.24884600000000001;  // deemax
    par[6] = 0.001;                // epsil
    par[7] = 1.149243;             // stmin
    auto medAir = new TGeoMedium("T_Air", numed, matAir, par);
    // Mixture: Other
    nel = 6;
    density = 3.2090000000000001;
    auto matOther = new TGeoMixture("Other", nel, density);
    a = 12.010999999999999;
    z = 6;
    w = 0.16870969999999999;  // C
    matOther->DefineElement(0, a, z, w);
    a = 1.0079400000000001;
    z = 1;
    w = 0.022652459999999999;  // H
    matOther->DefineElement(1, a, z, w);
    a = 15.999000000000001;
    z = 8;
    w = 0.089890460000000005;  // O
    matOther->DefineElement(2, a, z, w);
    a = 63.545999999999999;
    z = 29;
    w = 0.23839940000000001;  // CU
    matOther->DefineElement(3, a, z, w);
    a = 26.98;
    z = 13;
    w = 0.43784139999999999;  // AL
    matOther->DefineElement(4, a, z, w);
    a = 107.87;
    z = 47;
    w = 0.042506700000000001;  // AG
    matOther->DefineElement(5, a, z, w);
    matOther->SetIndex(88);
    // Medium: Other
    numed = 478;                    // medium number
    par[0] = 0;                     // isvol
    par[1] = -1;                    // ifield
    par[2] = 40;                    // fieldm
    par[3] = 20;                    // tmaxfd
    par[4] = 10000000000;           // stemax
    par[5] = 0.1710623;             // deemax
    par[6] = 0.5;                   // epsil
    par[7] = 0.032418339999999997;  // stmin
    auto medOther = new TGeoMedium("Pix_Bar_Cable", numed, matOther, par);

    // Material: T_Silicon
    a = 28.09;
    z = 14;
    density = 2.3300000000000001;
    radl = 9.3511055105325571;
    absl = 0;
    auto matSi = new TGeoMaterial("T_Silicon", a, z, density, radl, absl);
    matSi->SetIndex(76);

    // Medium: T_Silicon
    numed = 430;                   // medium number
    par[0] = 0;                    // isvol
    par[1] = -1;                   // ifield
    par[2] = 40;                   // fieldm
    par[3] = 20;                   // tmaxfd
    par[4] = 10000000000;          // stemax
    par[5] = 0.18459690000000001;  // deemax
    par[6] = 0.001;                // epsil
    par[7] = 0.03664369;           // stmin
    auto medSi = new TGeoMedium("T_Silicon", numed, matSi, par);

    // TRANSFORMATION MATRICES
    // Combi transformation:
    dx = 3.617388;
    dy = 2.0884999999999998;
    dz = 0;
    // Rotation: rot1306
    thx = 90;
    phx = 300;
    thy = 90;
    phy = 210;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix4 =
        new TGeoRotation("rot1306", thx, phx, thy, phy, thz, phz);
    auto pMatrix3 = new TGeoCombiTrans("", dx, dy, dz, pMatrix4);
    // Translation:
    dx = 0;
    dy = 0.058000000000000003;
    dz = 0;
    TGeoTranslation *pMatrix5 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0.0060000000000000001;
    dz = 16.68;
    TGeoTranslation *pMatrix6 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0.0060000000000000001;
    dz = -16.68;
    TGeoTranslation *pMatrix7 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0;
    dz = 20.010000000000002;
    TGeoTranslation *pMatrix8 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0;
    dz = -20.010000000000002;
    TGeoTranslation *pMatrix9 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = -0.0060000000000000001;
    dz = 23.34;
    TGeoTranslation *pMatrix10 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = -0.0060000000000000001;
    dz = -23.34;
    TGeoTranslation *pMatrix11 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = -0.0089999999999999993;
    dz = 0;
    TGeoTranslation *pMatrix12 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = -0.029000000000000001;
    dz = 0;
    TGeoTranslation *pMatrix14 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0.0089999999999999993;
    dz = 0;
    TGeoTranslation *pMatrix15 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0.033000000000000002;
    dz = 0;
    TGeoTranslation *pMatrix16 = new TGeoTranslation("", dx, dy, dz);
    // Translation:
    dx = 0;
    dy = 0.052999999999999999;
    dz = 0;
    TGeoTranslation *pMatrix17 = new TGeoTranslation("", dx, dy, dz);
    // Combi transformation:
    dx = 3.8269660000000001;
    dy = 2.2094999999999998;
    dz = 0;
    auto pMatrix18 = new TGeoCombiTrans("", dx, dy, dz, pMatrix4);
    // Combi transformation:
    dx = 2.975244;
    dy = 3.5457580000000002;
    dz = 0;
    // Rotation: rot1307
    thx = 90;
    phx = 140;
    thy = 90;
    phy = 50;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix20 =
        new TGeoRotation("rot1307", thx, phx, thy, phy, thz, phz);
    auto pMatrix19 = new TGeoCombiTrans("", dx, dy, dz, pMatrix20);
    // Combi transformation:
    dx = 2.81969;
    dy = 3.3603749999999999;
    dz = 0;
    auto pMatrix21 = new TGeoCombiTrans("", dx, dy, dz, pMatrix20);
    // Combi transformation:
    dx = 1.4286179999999999;
    dy = 3.9250959999999999;
    dz = 0;
    // Rotation: rot1308
    thx = 90;
    phx = 340;
    thy = 90;
    phy = 250;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix23 =
        new TGeoRotation("rot1308", thx, phx, thy, phy, thz, phz);
    auto pMatrix22 = new TGeoCombiTrans("", dx, dy, dz, pMatrix23);
    // Combi transformation:
    dx = 1.511387;
    dy = 4.1525020000000001;
    dz = 0;
    auto pMatrix24 = new TGeoCombiTrans("", dx, dy, dz, pMatrix23);
    // Combi transformation:
    dx = 5.6321349999999998e-08;
    dy = 4.6286579999999997;
    dz = 0;
    // Rotation: rot1309
    thx = 90;
    phx = 180;
    thy = 90;
    phy = 90;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix26 =
        new TGeoRotation("rot1309", thx, phx, thy, phy, thz, phz);
    auto pMatrix25 = new TGeoCombiTrans("", dx, dy, dz, pMatrix26);
    // Combi transformation:
    dx = 5.3376699999999999e-08;
    dy = 4.3866579999999997;
    dz = 0;
    auto pMatrix27 = new TGeoCombiTrans("", dx, dy, dz, pMatrix26);
    // Combi transformation:
    dx = -1.4286179999999999;
    dy = 3.9250959999999999;
    dz = 0;
    // Rotation: rot1310
    thx = 90;
    phx = 19.999999999999961;
    thy = 90;
    phy = 290;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix29 =
        new TGeoRotation("rot1310", thx, phx, thy, phy, thz, phz);
    auto pMatrix28 = new TGeoCombiTrans("", dx, dy, dz, pMatrix29);
    // Combi transformation:
    dx = -1.511387;
    dy = 4.1525020000000001;
    dz = 0;
    auto pMatrix30 = new TGeoCombiTrans("", dx, dy, dz, pMatrix29);
    // Combi transformation:
    dx = -2.975244;
    dy = 3.5457580000000002;
    dz = 0;
    // Rotation: rot1311
    thx = 90;
    phx = 219.99999999999997;
    thy = 90;
    phy = 130;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix32 =
        new TGeoRotation("rot1311", thx, phx, thy, phy, thz, phz);
    auto pMatrix31 = new TGeoCombiTrans("", dx, dy, dz, pMatrix32);
    // Combi transformation:
    dx = -2.8196889999999999;
    dy = 3.3603749999999999;
    dz = 0;
    auto pMatrix33 = new TGeoCombiTrans("", dx, dy, dz, pMatrix32);
    // Combi transformation:
    dx = -3.617388;
    dy = 2.0884999999999998;
    dz = 0;
    // Rotation: rot1312
    thx = 90;
    phx = 59.999999999999972;
    thy = 90;
    phy = 330;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix35 =
        new TGeoRotation("rot1312", thx, phx, thy, phy, thz, phz);
    auto pMatrix34 = new TGeoCombiTrans("", dx, dy, dz, pMatrix35);
    // Combi transformation:
    dx = -3.8269660000000001;
    dy = 2.2094999999999998;
    dz = 0;
    auto pMatrix36 = new TGeoCombiTrans("", dx, dy, dz, pMatrix35);
    // Combi transformation:
    dx = -4.5583390000000001;
    dy = 0.80375810000000003;
    dz = 0;
    // Rotation: rot1313
    thx = 90;
    phx = 260;
    thy = 90;
    phy = 170;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix38 =
        new TGeoRotation("rot1313", thx, phx, thy, phy, thz, phz);
    auto pMatrix37 = new TGeoCombiTrans("", dx, dy, dz, pMatrix38);
    // Combi transformation:
    dx = -4.3200149999999997;
    dy = 0.7617353;
    dz = 0;
    auto pMatrix39 = new TGeoCombiTrans("", dx, dy, dz, pMatrix38);
    // Combi transformation:
    dx = -4.1135419999999998;
    dy = -0.72532830000000004;
    dz = 0;
    // Rotation: rot1314
    thx = 90;
    phx = 99.999999999999972;
    thy = 90;
    phy = 9.9999999999999751;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix41 =
        new TGeoRotation("rot1314", thx, phx, thy, phy, thz, phz);
    auto pMatrix40 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);
    // Combi transformation:
    dx = -4.3518650000000001;
    dy = -0.76735120000000001;
    dz = 0;
    auto pMatrix42 = new TGeoCombiTrans("", dx, dy, dz, pMatrix41);
    // Combi transformation:
    dx = -4.0085350000000002;
    dy = -2.3143289999999999;
    dz = 0;
    // Rotation: rot1315
    thx = 90;
    phx = 300;
    thy = 90;
    phy = 210;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix44 =
        new TGeoRotation("rot1315", thx, phx, thy, phy, thz, phz);
    auto pMatrix43 = new TGeoCombiTrans("", dx, dy, dz, pMatrix44);
    // Combi transformation:
    dx = -3.7989570000000001;
    dy = -2.1933289999999999;
    dz = 0;
    auto pMatrix45 = new TGeoCombiTrans("", dx, dy, dz, pMatrix44);
    // Combi transformation:
    dx = -2.6849240000000001;
    dy = -3.1997680000000002;
    dz = 0;
    // Rotation: rot1316
    thx = 90;
    phx = 140.00000000000003;
    thy = 90;
    phy = 49.999999999999979;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix47 =
        new TGeoRotation("rot1316", thx, phx, thy, phy, thz, phz);
    auto pMatrix46 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);
    // Combi transformation:
    dx = -2.8404790000000002;
    dy = -3.3851499999999999;
    dz = 0;
    auto pMatrix48 = new TGeoCombiTrans("", dx, dy, dz, pMatrix47);
    // Combi transformation:
    dx = -1.5830949999999999;
    dy = -4.3495160000000004;
    dz = 0;
    // Rotation: rot1317
    thx = 90;
    phx = 340;
    thy = 90;
    phy = 250;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix50 =
        new TGeoRotation("rot1317", thx, phx, thy, phy, thz, phz);
    auto pMatrix49 = new TGeoCombiTrans("", dx, dy, dz, pMatrix50);
    // Combi transformation:
    dx = -1.500326;
    dy = -4.1221100000000002;
    dz = 0;
    auto pMatrix51 = new TGeoCombiTrans("", dx, dy, dz, pMatrix50);
    // Combi transformation:
    dx = -1.5247680000000001e-07;
    dy = -4.1769999999999996;
    dz = 0;
    // Rotation: rot1318
    thx = 90;
    phx = 180;
    thy = 90;
    phy = 90;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix53 =
        new TGeoRotation("rot1318", thx, phx, thy, phy, thz, phz);
    auto pMatrix52 = new TGeoCombiTrans("", dx, dy, dz, pMatrix53);
    // Combi transformation:
    dx = -1.613107e-07;
    dy = -4.4189999999999996;
    dz = 0;
    auto pMatrix54 = new TGeoCombiTrans("", dx, dy, dz, pMatrix53);
    // Combi transformation:
    dx = 1.583094;
    dy = -4.3495160000000004;
    dz = 0;
    // Rotation: rot1319
    thx = 90;
    phx = 19.999999999999961;
    thy = 90;
    phy = 290;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix56 =
        new TGeoRotation("rot1319", thx, phx, thy, phy, thz, phz);
    auto pMatrix55 = new TGeoCombiTrans("", dx, dy, dz, pMatrix56);
    // Combi transformation:
    dx = 1.5003249999999999;
    dy = -4.1221100000000002;
    dz = 0;
    auto pMatrix57 = new TGeoCombiTrans("", dx, dy, dz, pMatrix56);
    // Combi transformation:
    dx = 2.6849240000000001;
    dy = -3.1997680000000002;
    dz = 0;
    // Rotation: rot1320
    thx = 90;
    phx = 219.99999999999994;
    thy = 90;
    phy = 130.00000000000003;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix59 =
        new TGeoRotation("rot1320", thx, phx, thy, phy, thz, phz);
    auto pMatrix58 = new TGeoCombiTrans("", dx, dy, dz, pMatrix59);
    // Combi transformation:
    dx = 2.8404780000000001;
    dy = -3.385151;
    dz = 0;
    auto pMatrix60 = new TGeoCombiTrans("", dx, dy, dz, pMatrix59);
    // Combi transformation:
    dx = 4.0085350000000002;
    dy = -2.3143289999999999;
    dz = 0;
    // Rotation: rot1321
    thx = 90;
    phx = 59.999999999999972;
    thy = 90;
    phy = 330;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix62 =
        new TGeoRotation("rot1321", thx, phx, thy, phy, thz, phz);
    auto pMatrix61 = new TGeoCombiTrans("", dx, dy, dz, pMatrix62);
    // Combi transformation:
    dx = 3.7989570000000001;
    dy = -2.1933289999999999;
    dz = 0;
    auto pMatrix63 = new TGeoCombiTrans("", dx, dy, dz, pMatrix62);
    // Combi transformation:
    dx = 4.1135419999999998;
    dy = -0.72532859999999999;
    dz = 0;
    // Rotation: rot1322
    thx = 90;
    phx = 260;
    thy = 90;
    phy = 169.99999999999997;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix65 =
        new TGeoRotation("rot1322", thx, phx, thy, phy, thz, phz);
    auto pMatrix64 = new TGeoCombiTrans("", dx, dy, dz, pMatrix65);
    // Combi transformation:
    dx = 4.3518650000000001;
    dy = -0.76735149999999996;
    dz = 0;
    auto pMatrix66 = new TGeoCombiTrans("", dx, dy, dz, pMatrix65);
    // Combi transformation:
    dx = 4.5583390000000001;
    dy = 0.80375779999999997;
    dz = 0;
    // Rotation: rot1323
    thx = 90;
    phx = 99.999999999999972;
    thy = 90;
    phy = 9.9999999999999751;
    thz = 0;
    phz = 0;
    TGeoRotation *pMatrix68 =
        new TGeoRotation("rot1323", thx, phx, thy, phy, thz, phz);
    auto pMatrix67 = new TGeoCombiTrans("", dx, dy, dz, pMatrix68);
    // Combi transformation:
    dx = 4.3200149999999997;
    dy = 0.76173500000000005;
    dz = 0;
    auto pMatrix69 = new TGeoCombiTrans("", dx, dy, dz, pMatrix68);
    // Shape: TGeoPgon type: TGeoPgon
    phi1 = 0;
    dphi = 360;
    nedges = 18;
    nz = 2;
    TGeoPgon *pgon = new TGeoPgon("TGeoPgon", phi1, dphi, nedges, nz);
    z = -27.699999999999999;
    rmin = 4.1100000000000003;
    rmax = 4.8399999999999999;
    pgon->DefineSection(0, z, rmin, rmax);
    z = 27.699999999999999;
    rmin = 4.1100000000000003;
    rmax = 4.8399999999999999;
    pgon->DefineSection(1, z, rmin, rmax);
    TGeoShape *pTGeoPgon_682 = pgon;

    // Volume: PX1B
    auto pPX1B = new TGeoVolume(volumeName.c_str(), pTGeoPgon_682, medAir);

    // SET TOP VOLUME OF GEOMETRY
    gGeoManager->SetTopVolume(pPX1B);

    // SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.25;
    dy = 0.067000000000000004;
    dz = 26.68;
    TGeoShape *pTGeoBBox_681 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBL
    auto pPXBL = new TGeoVolume("PXBL", pTGeoBBox_681, medAir);
    pPX1B->AddNode(pPXBL, 1, pMatrix3);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 0.17499999999999999;
    dy = 0.17499999999999999;
    dz = 27.699999999999999;
    TGeoShape *pTGeoBBox_685 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBC

    auto pPXBC = new TGeoVolume("PXBC", pTGeoBBox_685, medOther);
    pPX1B->AddNode(pPXBL, 2, pMatrix19);
    pPX1B->AddNode(pPXBL, 3, pMatrix22);
    pPX1B->AddNode(pPXBL, 4, pMatrix25);
    pPX1B->AddNode(pPXBL, 5, pMatrix28);
    pPX1B->AddNode(pPXBL, 6, pMatrix31);
    pPX1B->AddNode(pPXBL, 7, pMatrix34);
    pPX1B->AddNode(pPXBL, 8, pMatrix37);
    pPX1B->AddNode(pPXBL, 9, pMatrix40);
    pPX1B->AddNode(pPXBL, 10, pMatrix43);
    pPX1B->AddNode(pPXBL, 11, pMatrix46);
    pPX1B->AddNode(pPXBL, 12, pMatrix49);
    pPX1B->AddNode(pPXBL, 13, pMatrix52);
    pPX1B->AddNode(pPXBL, 14, pMatrix55);
    pPX1B->AddNode(pPXBL, 15, pMatrix58);
    pPX1B->AddNode(pPXBL, 16, pMatrix61);
    pPX1B->AddNode(pPXBL, 17, pMatrix64);
    pPX1B->AddNode(pPXBL, 18, pMatrix67);

    pPX1B->AddNode(pPXBC, 1, pMatrix18);
    pPX1B->AddNode(pPXBC, 2, pMatrix21);
    pPX1B->AddNode(pPXBC, 3, pMatrix24);
    pPX1B->AddNode(pPXBC, 4, pMatrix27);
    pPX1B->AddNode(pPXBC, 5, pMatrix30);
    pPX1B->AddNode(pPXBC, 6, pMatrix33);
    pPX1B->AddNode(pPXBC, 7, pMatrix36);
    pPX1B->AddNode(pPXBC, 8, pMatrix39);
    pPX1B->AddNode(pPXBC, 9, pMatrix42);
    pPX1B->AddNode(pPXBC, 10, pMatrix45);
    pPX1B->AddNode(pPXBC, 11, pMatrix48);
    pPX1B->AddNode(pPXBC, 12, pMatrix51);
    pPX1B->AddNode(pPXBC, 13, pMatrix54);
    pPX1B->AddNode(pPXBC, 14, pMatrix57);
    pPX1B->AddNode(pPXBC, 15, pMatrix60);
    pPX1B->AddNode(pPXBC, 16, pMatrix63);
    pPX1B->AddNode(pPXBC, 17, pMatrix66);
    pPX1B->AddNode(pPXBC, 18, pMatrix69);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.25;
    dy = 0.0089999999999999993;
    dz = 26.68;
    TGeoShape *pTGeoBBox_679 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBB
    auto pPXBB_7fddef69dcd0 = new TGeoVolume("PXBB", pTGeoBBox_679, medAir);
    pPXBL->AddNode(pPXBB_7fddef69dcd0, 1, pMatrix5);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.25;
    dy = 0.058000000000000003;
    dz = 26.68;
    TGeoShape *pTGeoBBox_680 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBM
    auto pPXBM_7fddef7154b0 = new TGeoVolume("PXBM", pTGeoBBox_680, medAir);
    pPXBL->AddNode(pPXBM_7fddef7154b0, 1, pMatrix12);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.1000000000000001;
    dy = 0.0030000000000000001;
    dz = 10;
    TGeoShape *pTGeoBBox_686 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXB1
    auto pPXB1_7fddef69df10 = new TGeoVolume("PXB1", pTGeoBBox_686, medOther);
    pPXBB_7fddef69dcd0->AddNode(pPXB1_7fddef69df10, 1, pMatrix6);
    pPXBB_7fddef69dcd0->AddNode(pPXB1_7fddef69df10, 2, pMatrix7);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.1000000000000001;
    dy = 0.0030000000000000001;
    dz = 6.6699999999999999;
    TGeoShape *pTGeoBBox_687 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXB2
    auto pPXB2_7fddef69e460 = new TGeoVolume("PXB2", pTGeoBBox_687, medOther);
    pPXBB_7fddef69dcd0->AddNode(pPXB2_7fddef69e460, 1, pMatrix8);
    pPXBB_7fddef69dcd0->AddNode(pPXB2_7fddef69e460, 2, pMatrix9);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.1000000000000001;
    dy = 0.0030000000000000001;
    dz = 3.3399999999999999;
    TGeoShape *pTGeoBBox_688 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXB3
    auto pPXB3_7fddef69e670 = new TGeoVolume("PXB3", pTGeoBBox_688, medOther);
    pPXBB_7fddef69dcd0->AddNode(pPXB3_7fddef69e670, 1, pMatrix10);
    pPXBB_7fddef69dcd0->AddNode(pPXB3_7fddef69e670, 2, pMatrix11);
    TGeoVolume *pPXBQ_7fddef715760 =
        pPXBM_7fddef7154b0->Divide("PXBQ", 3, 8, -26.68, 6.6699999999999999);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.25;
    dy = 0.029000000000000001;
    dz = 3.3149999999999999;
    TGeoShape *pTGeoBBox_677 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBP
    auto pPXBP_7fddef715960 = new TGeoVolume("PXBP", pTGeoBBox_677, medOther);
    pPXBQ_7fddef715760->AddNode(pPXBP_7fddef715960, 1, pMatrix14);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 1.1000000000000001;
    dy = 0.0089999999999999993;
    dz = 3.2250000000000001;
    TGeoShape *pTGeoBBox_676 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBH
    auto pPXBH_7fddef715ed0 = new TGeoVolume("PXBH", pTGeoBBox_676, medOther);
    pPXBQ_7fddef715760->AddNode(pPXBH_7fddef715ed0, 1, pMatrix15);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 0.90000000000000002;
    dy = 0.014999999999999999;
    dz = 3.3149999999999999;
    TGeoShape *pTGeoBBox_674 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBE
    auto pPXBE_Silicon =
        new TGeoVolume(sensorName.c_str(), pTGeoBBox_674, medSi);
    pPXBQ_7fddef715760->AddNode(pPXBE_Silicon, 1, pMatrix16);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 0.81000000000000005;
    dy = 0.0050000000000000001;
    dz = 3.2000000000000002;
    TGeoShape *pTGeoBBox_678 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBA
    auto pPXBA_7fddef716d10 = new TGeoVolume("PXBA", pTGeoBBox_678, medOther);
    pPXBQ_7fddef715760->AddNode(pPXBA_7fddef716d10, 1, pMatrix17);
    // Shape: TGeoBBox type: TGeoBBox
    dx = 0.81000000000000005;
    dy = 0.014999999999999999;
    dz = 3.2250000000000001;
    TGeoShape *pTGeoBBox_675 = new TGeoBBox("TGeoBBox", dx, dy, dz);
    // Volume: PXBD
    auto pPXBD_7fddef716550 =
        new TGeoVolume(sensorName.c_str(), pTGeoBBox_675, medOther);
    pPXBE_Silicon->AddNode(pPXBD_7fddef716550, 1);

    // CLOSE GEOMETRY
    gGeoManager->CloseGeometry();
    return gGeoManager->GetTopVolume();
  }
};

}  // namespace Test

}  // namespace Acts