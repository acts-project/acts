//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<8, 8> multiply(const ActsMatrix<8, 1>& A, const ActsMatrix<1, 8>& B) {
  ActsMatrix<8, 8> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0];
  pOut[1] = pA[0]*pB[1];
  pOut[2] = pA[0]*pB[2];
  pOut[3] = pA[0]*pB[3];
  pOut[4] = pA[0]*pB[4];
  pOut[5] = pA[0]*pB[5];
  pOut[6] = pA[0]*pB[6];
  pOut[7] = pA[0]*pB[7];
  pOut[8] = pA[1]*pB[0];
  pOut[9] = pA[1]*pB[1];
  pOut[10] = pA[1]*pB[2];
  pOut[11] = pA[1]*pB[3];
  pOut[12] = pA[1]*pB[4];
  pOut[13] = pA[1]*pB[5];
  pOut[14] = pA[1]*pB[6];
  pOut[15] = pA[1]*pB[7];
  pOut[16] = pA[2]*pB[0];
  pOut[17] = pA[2]*pB[1];
  pOut[18] = pA[2]*pB[2];
  pOut[19] = pA[2]*pB[3];
  pOut[20] = pA[2]*pB[4];
  pOut[21] = pA[2]*pB[5];
  pOut[22] = pA[2]*pB[6];
  pOut[23] = pA[2]*pB[7];
  pOut[24] = pA[3]*pB[0];
  pOut[25] = pA[3]*pB[1];
  pOut[26] = pA[3]*pB[2];
  pOut[27] = pA[3]*pB[3];
  pOut[28] = pA[3]*pB[4];
  pOut[29] = pA[3]*pB[5];
  pOut[30] = pA[3]*pB[6];
  pOut[31] = pA[3]*pB[7];
  pOut[32] = pA[4]*pB[0];
  pOut[33] = pA[4]*pB[1];
  pOut[34] = pA[4]*pB[2];
  pOut[35] = pA[4]*pB[3];
  pOut[36] = pA[4]*pB[4];
  pOut[37] = pA[4]*pB[5];
  pOut[38] = pA[4]*pB[6];
  pOut[39] = pA[4]*pB[7];
  pOut[40] = pA[5]*pB[0];
  pOut[41] = pA[5]*pB[1];
  pOut[42] = pA[5]*pB[2];
  pOut[43] = pA[5]*pB[3];
  pOut[44] = pA[5]*pB[4];
  pOut[45] = pA[5]*pB[5];
  pOut[46] = pA[5]*pB[6];
  pOut[47] = pA[5]*pB[7];
  pOut[48] = pA[6]*pB[0];
  pOut[49] = pA[6]*pB[1];
  pOut[50] = pA[6]*pB[2];
  pOut[51] = pA[6]*pB[3];
  pOut[52] = pA[6]*pB[4];
  pOut[53] = pA[6]*pB[5];
  pOut[54] = pA[6]*pB[6];
  pOut[55] = pA[6]*pB[7];
  pOut[56] = pA[7]*pB[0];
  pOut[57] = pA[7]*pB[1];
  pOut[58] = pA[7]*pB[2];
  pOut[59] = pA[7]*pB[3];
  pOut[60] = pA[7]*pB[4];
  pOut[61] = pA[7]*pB[5];
  pOut[62] = pA[7]*pB[6];
  pOut[63] = pA[7]*pB[7];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<6, 8> multiply(const ActsMatrix<6, 8>& A, const ActsMatrix<8, 8>& B) {
  ActsMatrix<6, 8> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0] + pA[1]*pB[8] + pA[2]*pB[16] + pA[3]*pB[24] + pA[4]*pB[32] + pA[5]*pB[40] + pA[6]*pB[48] + pA[7]*pB[56];
  pOut[1] = pA[0]*pB[1] + pA[1]*pB[9] + pA[2]*pB[17] + pA[3]*pB[25] + pA[4]*pB[33] + pA[5]*pB[41] + pA[6]*pB[49] + pA[7]*pB[57];
  pOut[2] = pA[0]*pB[2] + pA[1]*pB[10] + pA[2]*pB[18] + pA[3]*pB[26] + pA[4]*pB[34] + pA[5]*pB[42] + pA[6]*pB[50] + pA[7]*pB[58];
  pOut[3] = pA[0]*pB[3] + pA[1]*pB[11] + pA[2]*pB[19] + pA[3]*pB[27] + pA[4]*pB[35] + pA[5]*pB[43] + pA[6]*pB[51] + pA[7]*pB[59];
  pOut[4] = pA[0]*pB[4] + pA[1]*pB[12] + pA[2]*pB[20] + pA[3]*pB[28] + pA[4]*pB[36] + pA[5]*pB[44] + pA[6]*pB[52] + pA[7]*pB[60];
  pOut[5] = pA[0]*pB[5] + pA[1]*pB[13] + pA[2]*pB[21] + pA[3]*pB[29] + pA[4]*pB[37] + pA[5]*pB[45] + pA[6]*pB[53] + pA[7]*pB[61];
  pOut[6] = pA[0]*pB[6] + pA[1]*pB[14] + pA[2]*pB[22] + pA[3]*pB[30] + pA[4]*pB[38] + pA[5]*pB[46] + pA[6]*pB[54] + pA[7]*pB[62];
  pOut[7] = pA[0]*pB[7] + pA[1]*pB[15] + pA[2]*pB[23] + pA[3]*pB[31] + pA[4]*pB[39] + pA[5]*pB[47] + pA[6]*pB[55] + pA[7]*pB[63];
  pOut[8] = pA[8]*pB[0] + pA[9]*pB[8] + pA[10]*pB[16] + pA[11]*pB[24] + pA[12]*pB[32] + pA[13]*pB[40] + pA[14]*pB[48] + pA[15]*pB[56];
  pOut[9] = pA[8]*pB[1] + pA[9]*pB[9] + pA[10]*pB[17] + pA[11]*pB[25] + pA[12]*pB[33] + pA[13]*pB[41] + pA[14]*pB[49] + pA[15]*pB[57];
  pOut[10] = pA[8]*pB[2] + pA[9]*pB[10] + pA[10]*pB[18] + pA[11]*pB[26] + pA[12]*pB[34] + pA[13]*pB[42] + pA[14]*pB[50] + pA[15]*pB[58];
  pOut[11] = pA[8]*pB[3] + pA[9]*pB[11] + pA[10]*pB[19] + pA[11]*pB[27] + pA[12]*pB[35] + pA[13]*pB[43] + pA[14]*pB[51] + pA[15]*pB[59];
  pOut[12] = pA[8]*pB[4] + pA[9]*pB[12] + pA[10]*pB[20] + pA[11]*pB[28] + pA[12]*pB[36] + pA[13]*pB[44] + pA[14]*pB[52] + pA[15]*pB[60];
  pOut[13] = pA[8]*pB[5] + pA[9]*pB[13] + pA[10]*pB[21] + pA[11]*pB[29] + pA[12]*pB[37] + pA[13]*pB[45] + pA[14]*pB[53] + pA[15]*pB[61];
  pOut[14] = pA[8]*pB[6] + pA[9]*pB[14] + pA[10]*pB[22] + pA[11]*pB[30] + pA[12]*pB[38] + pA[13]*pB[46] + pA[14]*pB[54] + pA[15]*pB[62];
  pOut[15] = pA[8]*pB[7] + pA[9]*pB[15] + pA[10]*pB[23] + pA[11]*pB[31] + pA[12]*pB[39] + pA[13]*pB[47] + pA[14]*pB[55] + pA[15]*pB[63];
  pOut[16] = pA[16]*pB[0] + pA[17]*pB[8] + pA[18]*pB[16] + pA[19]*pB[24] + pA[20]*pB[32] + pA[21]*pB[40] + pA[22]*pB[48] + pA[23]*pB[56];
  pOut[17] = pA[16]*pB[1] + pA[17]*pB[9] + pA[18]*pB[17] + pA[19]*pB[25] + pA[20]*pB[33] + pA[21]*pB[41] + pA[22]*pB[49] + pA[23]*pB[57];
  pOut[18] = pA[16]*pB[2] + pA[17]*pB[10] + pA[18]*pB[18] + pA[19]*pB[26] + pA[20]*pB[34] + pA[21]*pB[42] + pA[22]*pB[50] + pA[23]*pB[58];
  pOut[19] = pA[16]*pB[3] + pA[17]*pB[11] + pA[18]*pB[19] + pA[19]*pB[27] + pA[20]*pB[35] + pA[21]*pB[43] + pA[22]*pB[51] + pA[23]*pB[59];
  pOut[20] = pA[16]*pB[4] + pA[17]*pB[12] + pA[18]*pB[20] + pA[19]*pB[28] + pA[20]*pB[36] + pA[21]*pB[44] + pA[22]*pB[52] + pA[23]*pB[60];
  pOut[21] = pA[16]*pB[5] + pA[17]*pB[13] + pA[18]*pB[21] + pA[19]*pB[29] + pA[20]*pB[37] + pA[21]*pB[45] + pA[22]*pB[53] + pA[23]*pB[61];
  pOut[22] = pA[16]*pB[6] + pA[17]*pB[14] + pA[18]*pB[22] + pA[19]*pB[30] + pA[20]*pB[38] + pA[21]*pB[46] + pA[22]*pB[54] + pA[23]*pB[62];
  pOut[23] = pA[16]*pB[7] + pA[17]*pB[15] + pA[18]*pB[23] + pA[19]*pB[31] + pA[20]*pB[39] + pA[21]*pB[47] + pA[22]*pB[55] + pA[23]*pB[63];
  pOut[24] = pA[24]*pB[0] + pA[25]*pB[8] + pA[26]*pB[16] + pA[27]*pB[24] + pA[28]*pB[32] + pA[29]*pB[40] + pA[30]*pB[48] + pA[31]*pB[56];
  pOut[25] = pA[24]*pB[1] + pA[25]*pB[9] + pA[26]*pB[17] + pA[27]*pB[25] + pA[28]*pB[33] + pA[29]*pB[41] + pA[30]*pB[49] + pA[31]*pB[57];
  pOut[26] = pA[24]*pB[2] + pA[25]*pB[10] + pA[26]*pB[18] + pA[27]*pB[26] + pA[28]*pB[34] + pA[29]*pB[42] + pA[30]*pB[50] + pA[31]*pB[58];
  pOut[27] = pA[24]*pB[3] + pA[25]*pB[11] + pA[26]*pB[19] + pA[27]*pB[27] + pA[28]*pB[35] + pA[29]*pB[43] + pA[30]*pB[51] + pA[31]*pB[59];
  pOut[28] = pA[24]*pB[4] + pA[25]*pB[12] + pA[26]*pB[20] + pA[27]*pB[28] + pA[28]*pB[36] + pA[29]*pB[44] + pA[30]*pB[52] + pA[31]*pB[60];
  pOut[29] = pA[24]*pB[5] + pA[25]*pB[13] + pA[26]*pB[21] + pA[27]*pB[29] + pA[28]*pB[37] + pA[29]*pB[45] + pA[30]*pB[53] + pA[31]*pB[61];
  pOut[30] = pA[24]*pB[6] + pA[25]*pB[14] + pA[26]*pB[22] + pA[27]*pB[30] + pA[28]*pB[38] + pA[29]*pB[46] + pA[30]*pB[54] + pA[31]*pB[62];
  pOut[31] = pA[24]*pB[7] + pA[25]*pB[15] + pA[26]*pB[23] + pA[27]*pB[31] + pA[28]*pB[39] + pA[29]*pB[47] + pA[30]*pB[55] + pA[31]*pB[63];
  pOut[32] = pA[32]*pB[0] + pA[33]*pB[8] + pA[34]*pB[16] + pA[35]*pB[24] + pA[36]*pB[32] + pA[37]*pB[40] + pA[38]*pB[48] + pA[39]*pB[56];
  pOut[33] = pA[32]*pB[1] + pA[33]*pB[9] + pA[34]*pB[17] + pA[35]*pB[25] + pA[36]*pB[33] + pA[37]*pB[41] + pA[38]*pB[49] + pA[39]*pB[57];
  pOut[34] = pA[32]*pB[2] + pA[33]*pB[10] + pA[34]*pB[18] + pA[35]*pB[26] + pA[36]*pB[34] + pA[37]*pB[42] + pA[38]*pB[50] + pA[39]*pB[58];
  pOut[35] = pA[32]*pB[3] + pA[33]*pB[11] + pA[34]*pB[19] + pA[35]*pB[27] + pA[36]*pB[35] + pA[37]*pB[43] + pA[38]*pB[51] + pA[39]*pB[59];
  pOut[36] = pA[32]*pB[4] + pA[33]*pB[12] + pA[34]*pB[20] + pA[35]*pB[28] + pA[36]*pB[36] + pA[37]*pB[44] + pA[38]*pB[52] + pA[39]*pB[60];
  pOut[37] = pA[32]*pB[5] + pA[33]*pB[13] + pA[34]*pB[21] + pA[35]*pB[29] + pA[36]*pB[37] + pA[37]*pB[45] + pA[38]*pB[53] + pA[39]*pB[61];
  pOut[38] = pA[32]*pB[6] + pA[33]*pB[14] + pA[34]*pB[22] + pA[35]*pB[30] + pA[36]*pB[38] + pA[37]*pB[46] + pA[38]*pB[54] + pA[39]*pB[62];
  pOut[39] = pA[32]*pB[7] + pA[33]*pB[15] + pA[34]*pB[23] + pA[35]*pB[31] + pA[36]*pB[39] + pA[37]*pB[47] + pA[38]*pB[55] + pA[39]*pB[63];
  pOut[40] = pA[40]*pB[0] + pA[41]*pB[8] + pA[42]*pB[16] + pA[43]*pB[24] + pA[44]*pB[32] + pA[45]*pB[40] + pA[46]*pB[48] + pA[47]*pB[56];
  pOut[41] = pA[40]*pB[1] + pA[41]*pB[9] + pA[42]*pB[17] + pA[43]*pB[25] + pA[44]*pB[33] + pA[45]*pB[41] + pA[46]*pB[49] + pA[47]*pB[57];
  pOut[42] = pA[40]*pB[2] + pA[41]*pB[10] + pA[42]*pB[18] + pA[43]*pB[26] + pA[44]*pB[34] + pA[45]*pB[42] + pA[46]*pB[50] + pA[47]*pB[58];
  pOut[43] = pA[40]*pB[3] + pA[41]*pB[11] + pA[42]*pB[19] + pA[43]*pB[27] + pA[44]*pB[35] + pA[45]*pB[43] + pA[46]*pB[51] + pA[47]*pB[59];
  pOut[44] = pA[40]*pB[4] + pA[41]*pB[12] + pA[42]*pB[20] + pA[43]*pB[28] + pA[44]*pB[36] + pA[45]*pB[44] + pA[46]*pB[52] + pA[47]*pB[60];
  pOut[45] = pA[40]*pB[5] + pA[41]*pB[13] + pA[42]*pB[21] + pA[43]*pB[29] + pA[44]*pB[37] + pA[45]*pB[45] + pA[46]*pB[53] + pA[47]*pB[61];
  pOut[46] = pA[40]*pB[6] + pA[41]*pB[14] + pA[42]*pB[22] + pA[43]*pB[30] + pA[44]*pB[38] + pA[45]*pB[46] + pA[46]*pB[54] + pA[47]*pB[62];
  pOut[47] = pA[40]*pB[7] + pA[41]*pB[15] + pA[42]*pB[23] + pA[43]*pB[31] + pA[44]*pB[39] + pA[45]*pB[47] + pA[46]*pB[55] + pA[47]*pB[63];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<8, 8> multiply(const ActsMatrix<8, 8>& A, const ActsMatrix<8, 8>& B) {
  ActsMatrix<8, 8> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0] + pA[1]*pB[8] + pA[2]*pB[16] + pA[3]*pB[24] + pA[4]*pB[32] + pA[5]*pB[40] + pA[6]*pB[48] + pA[7]*pB[56];
  pOut[1] = pA[0]*pB[1] + pA[1]*pB[9] + pA[2]*pB[17] + pA[3]*pB[25] + pA[4]*pB[33] + pA[5]*pB[41] + pA[6]*pB[49] + pA[7]*pB[57];
  pOut[2] = pA[0]*pB[2] + pA[1]*pB[10] + pA[2]*pB[18] + pA[3]*pB[26] + pA[4]*pB[34] + pA[5]*pB[42] + pA[6]*pB[50] + pA[7]*pB[58];
  pOut[3] = pA[0]*pB[3] + pA[1]*pB[11] + pA[2]*pB[19] + pA[3]*pB[27] + pA[4]*pB[35] + pA[5]*pB[43] + pA[6]*pB[51] + pA[7]*pB[59];
  pOut[4] = pA[0]*pB[4] + pA[1]*pB[12] + pA[2]*pB[20] + pA[3]*pB[28] + pA[4]*pB[36] + pA[5]*pB[44] + pA[6]*pB[52] + pA[7]*pB[60];
  pOut[5] = pA[0]*pB[5] + pA[1]*pB[13] + pA[2]*pB[21] + pA[3]*pB[29] + pA[4]*pB[37] + pA[5]*pB[45] + pA[6]*pB[53] + pA[7]*pB[61];
  pOut[6] = pA[0]*pB[6] + pA[1]*pB[14] + pA[2]*pB[22] + pA[3]*pB[30] + pA[4]*pB[38] + pA[5]*pB[46] + pA[6]*pB[54] + pA[7]*pB[62];
  pOut[7] = pA[0]*pB[7] + pA[1]*pB[15] + pA[2]*pB[23] + pA[3]*pB[31] + pA[4]*pB[39] + pA[5]*pB[47] + pA[6]*pB[55] + pA[7]*pB[63];
  pOut[8] = pA[8]*pB[0] + pA[9]*pB[8] + pA[10]*pB[16] + pA[11]*pB[24] + pA[12]*pB[32] + pA[13]*pB[40] + pA[14]*pB[48] + pA[15]*pB[56];
  pOut[9] = pA[8]*pB[1] + pA[9]*pB[9] + pA[10]*pB[17] + pA[11]*pB[25] + pA[12]*pB[33] + pA[13]*pB[41] + pA[14]*pB[49] + pA[15]*pB[57];
  pOut[10] = pA[8]*pB[2] + pA[9]*pB[10] + pA[10]*pB[18] + pA[11]*pB[26] + pA[12]*pB[34] + pA[13]*pB[42] + pA[14]*pB[50] + pA[15]*pB[58];
  pOut[11] = pA[8]*pB[3] + pA[9]*pB[11] + pA[10]*pB[19] + pA[11]*pB[27] + pA[12]*pB[35] + pA[13]*pB[43] + pA[14]*pB[51] + pA[15]*pB[59];
  pOut[12] = pA[8]*pB[4] + pA[9]*pB[12] + pA[10]*pB[20] + pA[11]*pB[28] + pA[12]*pB[36] + pA[13]*pB[44] + pA[14]*pB[52] + pA[15]*pB[60];
  pOut[13] = pA[8]*pB[5] + pA[9]*pB[13] + pA[10]*pB[21] + pA[11]*pB[29] + pA[12]*pB[37] + pA[13]*pB[45] + pA[14]*pB[53] + pA[15]*pB[61];
  pOut[14] = pA[8]*pB[6] + pA[9]*pB[14] + pA[10]*pB[22] + pA[11]*pB[30] + pA[12]*pB[38] + pA[13]*pB[46] + pA[14]*pB[54] + pA[15]*pB[62];
  pOut[15] = pA[8]*pB[7] + pA[9]*pB[15] + pA[10]*pB[23] + pA[11]*pB[31] + pA[12]*pB[39] + pA[13]*pB[47] + pA[14]*pB[55] + pA[15]*pB[63];
  pOut[16] = pA[16]*pB[0] + pA[17]*pB[8] + pA[18]*pB[16] + pA[19]*pB[24] + pA[20]*pB[32] + pA[21]*pB[40] + pA[22]*pB[48] + pA[23]*pB[56];
  pOut[17] = pA[16]*pB[1] + pA[17]*pB[9] + pA[18]*pB[17] + pA[19]*pB[25] + pA[20]*pB[33] + pA[21]*pB[41] + pA[22]*pB[49] + pA[23]*pB[57];
  pOut[18] = pA[16]*pB[2] + pA[17]*pB[10] + pA[18]*pB[18] + pA[19]*pB[26] + pA[20]*pB[34] + pA[21]*pB[42] + pA[22]*pB[50] + pA[23]*pB[58];
  pOut[19] = pA[16]*pB[3] + pA[17]*pB[11] + pA[18]*pB[19] + pA[19]*pB[27] + pA[20]*pB[35] + pA[21]*pB[43] + pA[22]*pB[51] + pA[23]*pB[59];
  pOut[20] = pA[16]*pB[4] + pA[17]*pB[12] + pA[18]*pB[20] + pA[19]*pB[28] + pA[20]*pB[36] + pA[21]*pB[44] + pA[22]*pB[52] + pA[23]*pB[60];
  pOut[21] = pA[16]*pB[5] + pA[17]*pB[13] + pA[18]*pB[21] + pA[19]*pB[29] + pA[20]*pB[37] + pA[21]*pB[45] + pA[22]*pB[53] + pA[23]*pB[61];
  pOut[22] = pA[16]*pB[6] + pA[17]*pB[14] + pA[18]*pB[22] + pA[19]*pB[30] + pA[20]*pB[38] + pA[21]*pB[46] + pA[22]*pB[54] + pA[23]*pB[62];
  pOut[23] = pA[16]*pB[7] + pA[17]*pB[15] + pA[18]*pB[23] + pA[19]*pB[31] + pA[20]*pB[39] + pA[21]*pB[47] + pA[22]*pB[55] + pA[23]*pB[63];
  pOut[24] = pA[24]*pB[0] + pA[25]*pB[8] + pA[26]*pB[16] + pA[27]*pB[24] + pA[28]*pB[32] + pA[29]*pB[40] + pA[30]*pB[48] + pA[31]*pB[56];
  pOut[25] = pA[24]*pB[1] + pA[25]*pB[9] + pA[26]*pB[17] + pA[27]*pB[25] + pA[28]*pB[33] + pA[29]*pB[41] + pA[30]*pB[49] + pA[31]*pB[57];
  pOut[26] = pA[24]*pB[2] + pA[25]*pB[10] + pA[26]*pB[18] + pA[27]*pB[26] + pA[28]*pB[34] + pA[29]*pB[42] + pA[30]*pB[50] + pA[31]*pB[58];
  pOut[27] = pA[24]*pB[3] + pA[25]*pB[11] + pA[26]*pB[19] + pA[27]*pB[27] + pA[28]*pB[35] + pA[29]*pB[43] + pA[30]*pB[51] + pA[31]*pB[59];
  pOut[28] = pA[24]*pB[4] + pA[25]*pB[12] + pA[26]*pB[20] + pA[27]*pB[28] + pA[28]*pB[36] + pA[29]*pB[44] + pA[30]*pB[52] + pA[31]*pB[60];
  pOut[29] = pA[24]*pB[5] + pA[25]*pB[13] + pA[26]*pB[21] + pA[27]*pB[29] + pA[28]*pB[37] + pA[29]*pB[45] + pA[30]*pB[53] + pA[31]*pB[61];
  pOut[30] = pA[24]*pB[6] + pA[25]*pB[14] + pA[26]*pB[22] + pA[27]*pB[30] + pA[28]*pB[38] + pA[29]*pB[46] + pA[30]*pB[54] + pA[31]*pB[62];
  pOut[31] = pA[24]*pB[7] + pA[25]*pB[15] + pA[26]*pB[23] + pA[27]*pB[31] + pA[28]*pB[39] + pA[29]*pB[47] + pA[30]*pB[55] + pA[31]*pB[63];
  pOut[32] = pA[32]*pB[0] + pA[33]*pB[8] + pA[34]*pB[16] + pA[35]*pB[24] + pA[36]*pB[32] + pA[37]*pB[40] + pA[38]*pB[48] + pA[39]*pB[56];
  pOut[33] = pA[32]*pB[1] + pA[33]*pB[9] + pA[34]*pB[17] + pA[35]*pB[25] + pA[36]*pB[33] + pA[37]*pB[41] + pA[38]*pB[49] + pA[39]*pB[57];
  pOut[34] = pA[32]*pB[2] + pA[33]*pB[10] + pA[34]*pB[18] + pA[35]*pB[26] + pA[36]*pB[34] + pA[37]*pB[42] + pA[38]*pB[50] + pA[39]*pB[58];
  pOut[35] = pA[32]*pB[3] + pA[33]*pB[11] + pA[34]*pB[19] + pA[35]*pB[27] + pA[36]*pB[35] + pA[37]*pB[43] + pA[38]*pB[51] + pA[39]*pB[59];
  pOut[36] = pA[32]*pB[4] + pA[33]*pB[12] + pA[34]*pB[20] + pA[35]*pB[28] + pA[36]*pB[36] + pA[37]*pB[44] + pA[38]*pB[52] + pA[39]*pB[60];
  pOut[37] = pA[32]*pB[5] + pA[33]*pB[13] + pA[34]*pB[21] + pA[35]*pB[29] + pA[36]*pB[37] + pA[37]*pB[45] + pA[38]*pB[53] + pA[39]*pB[61];
  pOut[38] = pA[32]*pB[6] + pA[33]*pB[14] + pA[34]*pB[22] + pA[35]*pB[30] + pA[36]*pB[38] + pA[37]*pB[46] + pA[38]*pB[54] + pA[39]*pB[62];
  pOut[39] = pA[32]*pB[7] + pA[33]*pB[15] + pA[34]*pB[23] + pA[35]*pB[31] + pA[36]*pB[39] + pA[37]*pB[47] + pA[38]*pB[55] + pA[39]*pB[63];
  pOut[40] = pA[40]*pB[0] + pA[41]*pB[8] + pA[42]*pB[16] + pA[43]*pB[24] + pA[44]*pB[32] + pA[45]*pB[40] + pA[46]*pB[48] + pA[47]*pB[56];
  pOut[41] = pA[40]*pB[1] + pA[41]*pB[9] + pA[42]*pB[17] + pA[43]*pB[25] + pA[44]*pB[33] + pA[45]*pB[41] + pA[46]*pB[49] + pA[47]*pB[57];
  pOut[42] = pA[40]*pB[2] + pA[41]*pB[10] + pA[42]*pB[18] + pA[43]*pB[26] + pA[44]*pB[34] + pA[45]*pB[42] + pA[46]*pB[50] + pA[47]*pB[58];
  pOut[43] = pA[40]*pB[3] + pA[41]*pB[11] + pA[42]*pB[19] + pA[43]*pB[27] + pA[44]*pB[35] + pA[45]*pB[43] + pA[46]*pB[51] + pA[47]*pB[59];
  pOut[44] = pA[40]*pB[4] + pA[41]*pB[12] + pA[42]*pB[20] + pA[43]*pB[28] + pA[44]*pB[36] + pA[45]*pB[44] + pA[46]*pB[52] + pA[47]*pB[60];
  pOut[45] = pA[40]*pB[5] + pA[41]*pB[13] + pA[42]*pB[21] + pA[43]*pB[29] + pA[44]*pB[37] + pA[45]*pB[45] + pA[46]*pB[53] + pA[47]*pB[61];
  pOut[46] = pA[40]*pB[6] + pA[41]*pB[14] + pA[42]*pB[22] + pA[43]*pB[30] + pA[44]*pB[38] + pA[45]*pB[46] + pA[46]*pB[54] + pA[47]*pB[62];
  pOut[47] = pA[40]*pB[7] + pA[41]*pB[15] + pA[42]*pB[23] + pA[43]*pB[31] + pA[44]*pB[39] + pA[45]*pB[47] + pA[46]*pB[55] + pA[47]*pB[63];
  pOut[48] = pA[48]*pB[0] + pA[49]*pB[8] + pA[50]*pB[16] + pA[51]*pB[24] + pA[52]*pB[32] + pA[53]*pB[40] + pA[54]*pB[48] + pA[55]*pB[56];
  pOut[49] = pA[48]*pB[1] + pA[49]*pB[9] + pA[50]*pB[17] + pA[51]*pB[25] + pA[52]*pB[33] + pA[53]*pB[41] + pA[54]*pB[49] + pA[55]*pB[57];
  pOut[50] = pA[48]*pB[2] + pA[49]*pB[10] + pA[50]*pB[18] + pA[51]*pB[26] + pA[52]*pB[34] + pA[53]*pB[42] + pA[54]*pB[50] + pA[55]*pB[58];
  pOut[51] = pA[48]*pB[3] + pA[49]*pB[11] + pA[50]*pB[19] + pA[51]*pB[27] + pA[52]*pB[35] + pA[53]*pB[43] + pA[54]*pB[51] + pA[55]*pB[59];
  pOut[52] = pA[48]*pB[4] + pA[49]*pB[12] + pA[50]*pB[20] + pA[51]*pB[28] + pA[52]*pB[36] + pA[53]*pB[44] + pA[54]*pB[52] + pA[55]*pB[60];
  pOut[53] = pA[48]*pB[5] + pA[49]*pB[13] + pA[50]*pB[21] + pA[51]*pB[29] + pA[52]*pB[37] + pA[53]*pB[45] + pA[54]*pB[53] + pA[55]*pB[61];
  pOut[54] = pA[48]*pB[6] + pA[49]*pB[14] + pA[50]*pB[22] + pA[51]*pB[30] + pA[52]*pB[38] + pA[53]*pB[46] + pA[54]*pB[54] + pA[55]*pB[62];
  pOut[55] = pA[48]*pB[7] + pA[49]*pB[15] + pA[50]*pB[23] + pA[51]*pB[31] + pA[52]*pB[39] + pA[53]*pB[47] + pA[54]*pB[55] + pA[55]*pB[63];
  pOut[56] = pA[56]*pB[0] + pA[57]*pB[8] + pA[58]*pB[16] + pA[59]*pB[24] + pA[60]*pB[32] + pA[61]*pB[40] + pA[62]*pB[48] + pA[63]*pB[56];
  pOut[57] = pA[56]*pB[1] + pA[57]*pB[9] + pA[58]*pB[17] + pA[59]*pB[25] + pA[60]*pB[33] + pA[61]*pB[41] + pA[62]*pB[49] + pA[63]*pB[57];
  pOut[58] = pA[56]*pB[2] + pA[57]*pB[10] + pA[58]*pB[18] + pA[59]*pB[26] + pA[60]*pB[34] + pA[61]*pB[42] + pA[62]*pB[50] + pA[63]*pB[58];
  pOut[59] = pA[56]*pB[3] + pA[57]*pB[11] + pA[58]*pB[19] + pA[59]*pB[27] + pA[60]*pB[35] + pA[61]*pB[43] + pA[62]*pB[51] + pA[63]*pB[59];
  pOut[60] = pA[56]*pB[4] + pA[57]*pB[12] + pA[58]*pB[20] + pA[59]*pB[28] + pA[60]*pB[36] + pA[61]*pB[44] + pA[62]*pB[52] + pA[63]*pB[60];
  pOut[61] = pA[56]*pB[5] + pA[57]*pB[13] + pA[58]*pB[21] + pA[59]*pB[29] + pA[60]*pB[37] + pA[61]*pB[45] + pA[62]*pB[53] + pA[63]*pB[61];
  pOut[62] = pA[56]*pB[6] + pA[57]*pB[14] + pA[58]*pB[22] + pA[59]*pB[30] + pA[60]*pB[38] + pA[61]*pB[46] + pA[62]*pB[54] + pA[63]*pB[62];
  pOut[63] = pA[56]*pB[7] + pA[57]*pB[15] + pA[58]*pB[23] + pA[59]*pB[31] + pA[60]*pB[39] + pA[61]*pB[47] + pA[62]*pB[55] + pA[63]*pB[63];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<8, 6> multiply(const ActsMatrix<8, 8>& A, const ActsMatrix<8, 6>& B) {
  ActsMatrix<8, 6> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0] + pA[1]*pB[6] + pA[2]*pB[12] + pA[3]*pB[18] + pA[4]*pB[24] + pA[5]*pB[30] + pA[6]*pB[36] + pA[7]*pB[42];
  pOut[1] = pA[0]*pB[1] + pA[1]*pB[7] + pA[2]*pB[13] + pA[3]*pB[19] + pA[4]*pB[25] + pA[5]*pB[31] + pA[6]*pB[37] + pA[7]*pB[43];
  pOut[2] = pA[0]*pB[2] + pA[1]*pB[8] + pA[2]*pB[14] + pA[3]*pB[20] + pA[4]*pB[26] + pA[5]*pB[32] + pA[6]*pB[38] + pA[7]*pB[44];
  pOut[3] = pA[0]*pB[3] + pA[1]*pB[9] + pA[2]*pB[15] + pA[3]*pB[21] + pA[4]*pB[27] + pA[5]*pB[33] + pA[6]*pB[39] + pA[7]*pB[45];
  pOut[4] = pA[0]*pB[4] + pA[1]*pB[10] + pA[2]*pB[16] + pA[3]*pB[22] + pA[4]*pB[28] + pA[5]*pB[34] + pA[6]*pB[40] + pA[7]*pB[46];
  pOut[5] = pA[0]*pB[5] + pA[1]*pB[11] + pA[2]*pB[17] + pA[3]*pB[23] + pA[4]*pB[29] + pA[5]*pB[35] + pA[6]*pB[41] + pA[7]*pB[47];
  pOut[6] = pA[8]*pB[0] + pA[9]*pB[6] + pA[10]*pB[12] + pA[11]*pB[18] + pA[12]*pB[24] + pA[13]*pB[30] + pA[14]*pB[36] + pA[15]*pB[42];
  pOut[7] = pA[8]*pB[1] + pA[9]*pB[7] + pA[10]*pB[13] + pA[11]*pB[19] + pA[12]*pB[25] + pA[13]*pB[31] + pA[14]*pB[37] + pA[15]*pB[43];
  pOut[8] = pA[8]*pB[2] + pA[9]*pB[8] + pA[10]*pB[14] + pA[11]*pB[20] + pA[12]*pB[26] + pA[13]*pB[32] + pA[14]*pB[38] + pA[15]*pB[44];
  pOut[9] = pA[8]*pB[3] + pA[9]*pB[9] + pA[10]*pB[15] + pA[11]*pB[21] + pA[12]*pB[27] + pA[13]*pB[33] + pA[14]*pB[39] + pA[15]*pB[45];
  pOut[10] = pA[8]*pB[4] + pA[9]*pB[10] + pA[10]*pB[16] + pA[11]*pB[22] + pA[12]*pB[28] + pA[13]*pB[34] + pA[14]*pB[40] + pA[15]*pB[46];
  pOut[11] = pA[8]*pB[5] + pA[9]*pB[11] + pA[10]*pB[17] + pA[11]*pB[23] + pA[12]*pB[29] + pA[13]*pB[35] + pA[14]*pB[41] + pA[15]*pB[47];
  pOut[12] = pA[16]*pB[0] + pA[17]*pB[6] + pA[18]*pB[12] + pA[19]*pB[18] + pA[20]*pB[24] + pA[21]*pB[30] + pA[22]*pB[36] + pA[23]*pB[42];
  pOut[13] = pA[16]*pB[1] + pA[17]*pB[7] + pA[18]*pB[13] + pA[19]*pB[19] + pA[20]*pB[25] + pA[21]*pB[31] + pA[22]*pB[37] + pA[23]*pB[43];
  pOut[14] = pA[16]*pB[2] + pA[17]*pB[8] + pA[18]*pB[14] + pA[19]*pB[20] + pA[20]*pB[26] + pA[21]*pB[32] + pA[22]*pB[38] + pA[23]*pB[44];
  pOut[15] = pA[16]*pB[3] + pA[17]*pB[9] + pA[18]*pB[15] + pA[19]*pB[21] + pA[20]*pB[27] + pA[21]*pB[33] + pA[22]*pB[39] + pA[23]*pB[45];
  pOut[16] = pA[16]*pB[4] + pA[17]*pB[10] + pA[18]*pB[16] + pA[19]*pB[22] + pA[20]*pB[28] + pA[21]*pB[34] + pA[22]*pB[40] + pA[23]*pB[46];
  pOut[17] = pA[16]*pB[5] + pA[17]*pB[11] + pA[18]*pB[17] + pA[19]*pB[23] + pA[20]*pB[29] + pA[21]*pB[35] + pA[22]*pB[41] + pA[23]*pB[47];
  pOut[18] = pA[24]*pB[0] + pA[25]*pB[6] + pA[26]*pB[12] + pA[27]*pB[18] + pA[28]*pB[24] + pA[29]*pB[30] + pA[30]*pB[36] + pA[31]*pB[42];
  pOut[19] = pA[24]*pB[1] + pA[25]*pB[7] + pA[26]*pB[13] + pA[27]*pB[19] + pA[28]*pB[25] + pA[29]*pB[31] + pA[30]*pB[37] + pA[31]*pB[43];
  pOut[20] = pA[24]*pB[2] + pA[25]*pB[8] + pA[26]*pB[14] + pA[27]*pB[20] + pA[28]*pB[26] + pA[29]*pB[32] + pA[30]*pB[38] + pA[31]*pB[44];
  pOut[21] = pA[24]*pB[3] + pA[25]*pB[9] + pA[26]*pB[15] + pA[27]*pB[21] + pA[28]*pB[27] + pA[29]*pB[33] + pA[30]*pB[39] + pA[31]*pB[45];
  pOut[22] = pA[24]*pB[4] + pA[25]*pB[10] + pA[26]*pB[16] + pA[27]*pB[22] + pA[28]*pB[28] + pA[29]*pB[34] + pA[30]*pB[40] + pA[31]*pB[46];
  pOut[23] = pA[24]*pB[5] + pA[25]*pB[11] + pA[26]*pB[17] + pA[27]*pB[23] + pA[28]*pB[29] + pA[29]*pB[35] + pA[30]*pB[41] + pA[31]*pB[47];
  pOut[24] = pA[32]*pB[0] + pA[33]*pB[6] + pA[34]*pB[12] + pA[35]*pB[18] + pA[36]*pB[24] + pA[37]*pB[30] + pA[38]*pB[36] + pA[39]*pB[42];
  pOut[25] = pA[32]*pB[1] + pA[33]*pB[7] + pA[34]*pB[13] + pA[35]*pB[19] + pA[36]*pB[25] + pA[37]*pB[31] + pA[38]*pB[37] + pA[39]*pB[43];
  pOut[26] = pA[32]*pB[2] + pA[33]*pB[8] + pA[34]*pB[14] + pA[35]*pB[20] + pA[36]*pB[26] + pA[37]*pB[32] + pA[38]*pB[38] + pA[39]*pB[44];
  pOut[27] = pA[32]*pB[3] + pA[33]*pB[9] + pA[34]*pB[15] + pA[35]*pB[21] + pA[36]*pB[27] + pA[37]*pB[33] + pA[38]*pB[39] + pA[39]*pB[45];
  pOut[28] = pA[32]*pB[4] + pA[33]*pB[10] + pA[34]*pB[16] + pA[35]*pB[22] + pA[36]*pB[28] + pA[37]*pB[34] + pA[38]*pB[40] + pA[39]*pB[46];
  pOut[29] = pA[32]*pB[5] + pA[33]*pB[11] + pA[34]*pB[17] + pA[35]*pB[23] + pA[36]*pB[29] + pA[37]*pB[35] + pA[38]*pB[41] + pA[39]*pB[47];
  pOut[30] = pA[40]*pB[0] + pA[41]*pB[6] + pA[42]*pB[12] + pA[43]*pB[18] + pA[44]*pB[24] + pA[45]*pB[30] + pA[46]*pB[36] + pA[47]*pB[42];
  pOut[31] = pA[40]*pB[1] + pA[41]*pB[7] + pA[42]*pB[13] + pA[43]*pB[19] + pA[44]*pB[25] + pA[45]*pB[31] + pA[46]*pB[37] + pA[47]*pB[43];
  pOut[32] = pA[40]*pB[2] + pA[41]*pB[8] + pA[42]*pB[14] + pA[43]*pB[20] + pA[44]*pB[26] + pA[45]*pB[32] + pA[46]*pB[38] + pA[47]*pB[44];
  pOut[33] = pA[40]*pB[3] + pA[41]*pB[9] + pA[42]*pB[15] + pA[43]*pB[21] + pA[44]*pB[27] + pA[45]*pB[33] + pA[46]*pB[39] + pA[47]*pB[45];
  pOut[34] = pA[40]*pB[4] + pA[41]*pB[10] + pA[42]*pB[16] + pA[43]*pB[22] + pA[44]*pB[28] + pA[45]*pB[34] + pA[46]*pB[40] + pA[47]*pB[46];
  pOut[35] = pA[40]*pB[5] + pA[41]*pB[11] + pA[42]*pB[17] + pA[43]*pB[23] + pA[44]*pB[29] + pA[45]*pB[35] + pA[46]*pB[41] + pA[47]*pB[47];
  pOut[36] = pA[48]*pB[0] + pA[49]*pB[6] + pA[50]*pB[12] + pA[51]*pB[18] + pA[52]*pB[24] + pA[53]*pB[30] + pA[54]*pB[36] + pA[55]*pB[42];
  pOut[37] = pA[48]*pB[1] + pA[49]*pB[7] + pA[50]*pB[13] + pA[51]*pB[19] + pA[52]*pB[25] + pA[53]*pB[31] + pA[54]*pB[37] + pA[55]*pB[43];
  pOut[38] = pA[48]*pB[2] + pA[49]*pB[8] + pA[50]*pB[14] + pA[51]*pB[20] + pA[52]*pB[26] + pA[53]*pB[32] + pA[54]*pB[38] + pA[55]*pB[44];
  pOut[39] = pA[48]*pB[3] + pA[49]*pB[9] + pA[50]*pB[15] + pA[51]*pB[21] + pA[52]*pB[27] + pA[53]*pB[33] + pA[54]*pB[39] + pA[55]*pB[45];
  pOut[40] = pA[48]*pB[4] + pA[49]*pB[10] + pA[50]*pB[16] + pA[51]*pB[22] + pA[52]*pB[28] + pA[53]*pB[34] + pA[54]*pB[40] + pA[55]*pB[46];
  pOut[41] = pA[48]*pB[5] + pA[49]*pB[11] + pA[50]*pB[17] + pA[51]*pB[23] + pA[52]*pB[29] + pA[53]*pB[35] + pA[54]*pB[41] + pA[55]*pB[47];
  pOut[42] = pA[56]*pB[0] + pA[57]*pB[6] + pA[58]*pB[12] + pA[59]*pB[18] + pA[60]*pB[24] + pA[61]*pB[30] + pA[62]*pB[36] + pA[63]*pB[42];
  pOut[43] = pA[56]*pB[1] + pA[57]*pB[7] + pA[58]*pB[13] + pA[59]*pB[19] + pA[60]*pB[25] + pA[61]*pB[31] + pA[62]*pB[37] + pA[63]*pB[43];
  pOut[44] = pA[56]*pB[2] + pA[57]*pB[8] + pA[58]*pB[14] + pA[59]*pB[20] + pA[60]*pB[26] + pA[61]*pB[32] + pA[62]*pB[38] + pA[63]*pB[44];
  pOut[45] = pA[56]*pB[3] + pA[57]*pB[9] + pA[58]*pB[15] + pA[59]*pB[21] + pA[60]*pB[27] + pA[61]*pB[33] + pA[62]*pB[39] + pA[63]*pB[45];
  pOut[46] = pA[56]*pB[4] + pA[57]*pB[10] + pA[58]*pB[16] + pA[59]*pB[22] + pA[60]*pB[28] + pA[61]*pB[34] + pA[62]*pB[40] + pA[63]*pB[46];
  pOut[47] = pA[56]*pB[5] + pA[57]*pB[11] + pA[58]*pB[17] + pA[59]*pB[23] + pA[60]*pB[29] + pA[61]*pB[35] + pA[62]*pB[41] + pA[63]*pB[47];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<6, 6> multiply(const ActsMatrix<6, 8>& A, const ActsMatrix<8, 6>& B) {
  ActsMatrix<6, 6> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0] + pA[1]*pB[6] + pA[2]*pB[12] + pA[3]*pB[18] + pA[4]*pB[24] + pA[5]*pB[30] + pA[6]*pB[36] + pA[7]*pB[42];
  pOut[1] = pA[0]*pB[1] + pA[1]*pB[7] + pA[2]*pB[13] + pA[3]*pB[19] + pA[4]*pB[25] + pA[5]*pB[31] + pA[6]*pB[37] + pA[7]*pB[43];
  pOut[2] = pA[0]*pB[2] + pA[1]*pB[8] + pA[2]*pB[14] + pA[3]*pB[20] + pA[4]*pB[26] + pA[5]*pB[32] + pA[6]*pB[38] + pA[7]*pB[44];
  pOut[3] = pA[0]*pB[3] + pA[1]*pB[9] + pA[2]*pB[15] + pA[3]*pB[21] + pA[4]*pB[27] + pA[5]*pB[33] + pA[6]*pB[39] + pA[7]*pB[45];
  pOut[4] = pA[0]*pB[4] + pA[1]*pB[10] + pA[2]*pB[16] + pA[3]*pB[22] + pA[4]*pB[28] + pA[5]*pB[34] + pA[6]*pB[40] + pA[7]*pB[46];
  pOut[5] = pA[0]*pB[5] + pA[1]*pB[11] + pA[2]*pB[17] + pA[3]*pB[23] + pA[4]*pB[29] + pA[5]*pB[35] + pA[6]*pB[41] + pA[7]*pB[47];
  pOut[6] = pA[8]*pB[0] + pA[9]*pB[6] + pA[10]*pB[12] + pA[11]*pB[18] + pA[12]*pB[24] + pA[13]*pB[30] + pA[14]*pB[36] + pA[15]*pB[42];
  pOut[7] = pA[8]*pB[1] + pA[9]*pB[7] + pA[10]*pB[13] + pA[11]*pB[19] + pA[12]*pB[25] + pA[13]*pB[31] + pA[14]*pB[37] + pA[15]*pB[43];
  pOut[8] = pA[8]*pB[2] + pA[9]*pB[8] + pA[10]*pB[14] + pA[11]*pB[20] + pA[12]*pB[26] + pA[13]*pB[32] + pA[14]*pB[38] + pA[15]*pB[44];
  pOut[9] = pA[8]*pB[3] + pA[9]*pB[9] + pA[10]*pB[15] + pA[11]*pB[21] + pA[12]*pB[27] + pA[13]*pB[33] + pA[14]*pB[39] + pA[15]*pB[45];
  pOut[10] = pA[8]*pB[4] + pA[9]*pB[10] + pA[10]*pB[16] + pA[11]*pB[22] + pA[12]*pB[28] + pA[13]*pB[34] + pA[14]*pB[40] + pA[15]*pB[46];
  pOut[11] = pA[8]*pB[5] + pA[9]*pB[11] + pA[10]*pB[17] + pA[11]*pB[23] + pA[12]*pB[29] + pA[13]*pB[35] + pA[14]*pB[41] + pA[15]*pB[47];
  pOut[12] = pA[16]*pB[0] + pA[17]*pB[6] + pA[18]*pB[12] + pA[19]*pB[18] + pA[20]*pB[24] + pA[21]*pB[30] + pA[22]*pB[36] + pA[23]*pB[42];
  pOut[13] = pA[16]*pB[1] + pA[17]*pB[7] + pA[18]*pB[13] + pA[19]*pB[19] + pA[20]*pB[25] + pA[21]*pB[31] + pA[22]*pB[37] + pA[23]*pB[43];
  pOut[14] = pA[16]*pB[2] + pA[17]*pB[8] + pA[18]*pB[14] + pA[19]*pB[20] + pA[20]*pB[26] + pA[21]*pB[32] + pA[22]*pB[38] + pA[23]*pB[44];
  pOut[15] = pA[16]*pB[3] + pA[17]*pB[9] + pA[18]*pB[15] + pA[19]*pB[21] + pA[20]*pB[27] + pA[21]*pB[33] + pA[22]*pB[39] + pA[23]*pB[45];
  pOut[16] = pA[16]*pB[4] + pA[17]*pB[10] + pA[18]*pB[16] + pA[19]*pB[22] + pA[20]*pB[28] + pA[21]*pB[34] + pA[22]*pB[40] + pA[23]*pB[46];
  pOut[17] = pA[16]*pB[5] + pA[17]*pB[11] + pA[18]*pB[17] + pA[19]*pB[23] + pA[20]*pB[29] + pA[21]*pB[35] + pA[22]*pB[41] + pA[23]*pB[47];
  pOut[18] = pA[24]*pB[0] + pA[25]*pB[6] + pA[26]*pB[12] + pA[27]*pB[18] + pA[28]*pB[24] + pA[29]*pB[30] + pA[30]*pB[36] + pA[31]*pB[42];
  pOut[19] = pA[24]*pB[1] + pA[25]*pB[7] + pA[26]*pB[13] + pA[27]*pB[19] + pA[28]*pB[25] + pA[29]*pB[31] + pA[30]*pB[37] + pA[31]*pB[43];
  pOut[20] = pA[24]*pB[2] + pA[25]*pB[8] + pA[26]*pB[14] + pA[27]*pB[20] + pA[28]*pB[26] + pA[29]*pB[32] + pA[30]*pB[38] + pA[31]*pB[44];
  pOut[21] = pA[24]*pB[3] + pA[25]*pB[9] + pA[26]*pB[15] + pA[27]*pB[21] + pA[28]*pB[27] + pA[29]*pB[33] + pA[30]*pB[39] + pA[31]*pB[45];
  pOut[22] = pA[24]*pB[4] + pA[25]*pB[10] + pA[26]*pB[16] + pA[27]*pB[22] + pA[28]*pB[28] + pA[29]*pB[34] + pA[30]*pB[40] + pA[31]*pB[46];
  pOut[23] = pA[24]*pB[5] + pA[25]*pB[11] + pA[26]*pB[17] + pA[27]*pB[23] + pA[28]*pB[29] + pA[29]*pB[35] + pA[30]*pB[41] + pA[31]*pB[47];
  pOut[24] = pA[32]*pB[0] + pA[33]*pB[6] + pA[34]*pB[12] + pA[35]*pB[18] + pA[36]*pB[24] + pA[37]*pB[30] + pA[38]*pB[36] + pA[39]*pB[42];
  pOut[25] = pA[32]*pB[1] + pA[33]*pB[7] + pA[34]*pB[13] + pA[35]*pB[19] + pA[36]*pB[25] + pA[37]*pB[31] + pA[38]*pB[37] + pA[39]*pB[43];
  pOut[26] = pA[32]*pB[2] + pA[33]*pB[8] + pA[34]*pB[14] + pA[35]*pB[20] + pA[36]*pB[26] + pA[37]*pB[32] + pA[38]*pB[38] + pA[39]*pB[44];
  pOut[27] = pA[32]*pB[3] + pA[33]*pB[9] + pA[34]*pB[15] + pA[35]*pB[21] + pA[36]*pB[27] + pA[37]*pB[33] + pA[38]*pB[39] + pA[39]*pB[45];
  pOut[28] = pA[32]*pB[4] + pA[33]*pB[10] + pA[34]*pB[16] + pA[35]*pB[22] + pA[36]*pB[28] + pA[37]*pB[34] + pA[38]*pB[40] + pA[39]*pB[46];
  pOut[29] = pA[32]*pB[5] + pA[33]*pB[11] + pA[34]*pB[17] + pA[35]*pB[23] + pA[36]*pB[29] + pA[37]*pB[35] + pA[38]*pB[41] + pA[39]*pB[47];
  pOut[30] = pA[40]*pB[0] + pA[41]*pB[6] + pA[42]*pB[12] + pA[43]*pB[18] + pA[44]*pB[24] + pA[45]*pB[30] + pA[46]*pB[36] + pA[47]*pB[42];
  pOut[31] = pA[40]*pB[1] + pA[41]*pB[7] + pA[42]*pB[13] + pA[43]*pB[19] + pA[44]*pB[25] + pA[45]*pB[31] + pA[46]*pB[37] + pA[47]*pB[43];
  pOut[32] = pA[40]*pB[2] + pA[41]*pB[8] + pA[42]*pB[14] + pA[43]*pB[20] + pA[44]*pB[26] + pA[45]*pB[32] + pA[46]*pB[38] + pA[47]*pB[44];
  pOut[33] = pA[40]*pB[3] + pA[41]*pB[9] + pA[42]*pB[15] + pA[43]*pB[21] + pA[44]*pB[27] + pA[45]*pB[33] + pA[46]*pB[39] + pA[47]*pB[45];
  pOut[34] = pA[40]*pB[4] + pA[41]*pB[10] + pA[42]*pB[16] + pA[43]*pB[22] + pA[44]*pB[28] + pA[45]*pB[34] + pA[46]*pB[40] + pA[47]*pB[46];
  pOut[35] = pA[40]*pB[5] + pA[41]*pB[11] + pA[42]*pB[17] + pA[43]*pB[23] + pA[44]*pB[29] + pA[45]*pB[35] + pA[46]*pB[41] + pA[47]*pB[47];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<6, 6> multiply(const ActsMatrix<6, 6>& A, const ActsMatrix<6, 6>& B) {
  ActsMatrix<6, 6> out;
  double* pOut = out.data();

  const double* pA = A.data();
  const double* pB = B.data();

  pOut[0] = pA[0]*pB[0] + pA[1]*pB[6] + pA[2]*pB[12] + pA[3]*pB[18] + pA[4]*pB[24] + pA[5]*pB[30];
  pOut[1] = pA[0]*pB[1] + pA[1]*pB[7] + pA[2]*pB[13] + pA[3]*pB[19] + pA[4]*pB[25] + pA[5]*pB[31];
  pOut[2] = pA[0]*pB[2] + pA[1]*pB[8] + pA[2]*pB[14] + pA[3]*pB[20] + pA[4]*pB[26] + pA[5]*pB[32];
  pOut[3] = pA[0]*pB[3] + pA[1]*pB[9] + pA[2]*pB[15] + pA[3]*pB[21] + pA[4]*pB[27] + pA[5]*pB[33];
  pOut[4] = pA[0]*pB[4] + pA[1]*pB[10] + pA[2]*pB[16] + pA[3]*pB[22] + pA[4]*pB[28] + pA[5]*pB[34];
  pOut[5] = pA[0]*pB[5] + pA[1]*pB[11] + pA[2]*pB[17] + pA[3]*pB[23] + pA[4]*pB[29] + pA[5]*pB[35];
  pOut[6] = pA[6]*pB[0] + pA[7]*pB[6] + pA[8]*pB[12] + pA[9]*pB[18] + pA[10]*pB[24] + pA[11]*pB[30];
  pOut[7] = pA[6]*pB[1] + pA[7]*pB[7] + pA[8]*pB[13] + pA[9]*pB[19] + pA[10]*pB[25] + pA[11]*pB[31];
  pOut[8] = pA[6]*pB[2] + pA[7]*pB[8] + pA[8]*pB[14] + pA[9]*pB[20] + pA[10]*pB[26] + pA[11]*pB[32];
  pOut[9] = pA[6]*pB[3] + pA[7]*pB[9] + pA[8]*pB[15] + pA[9]*pB[21] + pA[10]*pB[27] + pA[11]*pB[33];
  pOut[10] = pA[6]*pB[4] + pA[7]*pB[10] + pA[8]*pB[16] + pA[9]*pB[22] + pA[10]*pB[28] + pA[11]*pB[34];
  pOut[11] = pA[6]*pB[5] + pA[7]*pB[11] + pA[8]*pB[17] + pA[9]*pB[23] + pA[10]*pB[29] + pA[11]*pB[35];
  pOut[12] = pA[12]*pB[0] + pA[13]*pB[6] + pA[14]*pB[12] + pA[15]*pB[18] + pA[16]*pB[24] + pA[17]*pB[30];
  pOut[13] = pA[12]*pB[1] + pA[13]*pB[7] + pA[14]*pB[13] + pA[15]*pB[19] + pA[16]*pB[25] + pA[17]*pB[31];
  pOut[14] = pA[12]*pB[2] + pA[13]*pB[8] + pA[14]*pB[14] + pA[15]*pB[20] + pA[16]*pB[26] + pA[17]*pB[32];
  pOut[15] = pA[12]*pB[3] + pA[13]*pB[9] + pA[14]*pB[15] + pA[15]*pB[21] + pA[16]*pB[27] + pA[17]*pB[33];
  pOut[16] = pA[12]*pB[4] + pA[13]*pB[10] + pA[14]*pB[16] + pA[15]*pB[22] + pA[16]*pB[28] + pA[17]*pB[34];
  pOut[17] = pA[12]*pB[5] + pA[13]*pB[11] + pA[14]*pB[17] + pA[15]*pB[23] + pA[16]*pB[29] + pA[17]*pB[35];
  pOut[18] = pA[18]*pB[0] + pA[19]*pB[6] + pA[20]*pB[12] + pA[21]*pB[18] + pA[22]*pB[24] + pA[23]*pB[30];
  pOut[19] = pA[18]*pB[1] + pA[19]*pB[7] + pA[20]*pB[13] + pA[21]*pB[19] + pA[22]*pB[25] + pA[23]*pB[31];
  pOut[20] = pA[18]*pB[2] + pA[19]*pB[8] + pA[20]*pB[14] + pA[21]*pB[20] + pA[22]*pB[26] + pA[23]*pB[32];
  pOut[21] = pA[18]*pB[3] + pA[19]*pB[9] + pA[20]*pB[15] + pA[21]*pB[21] + pA[22]*pB[27] + pA[23]*pB[33];
  pOut[22] = pA[18]*pB[4] + pA[19]*pB[10] + pA[20]*pB[16] + pA[21]*pB[22] + pA[22]*pB[28] + pA[23]*pB[34];
  pOut[23] = pA[18]*pB[5] + pA[19]*pB[11] + pA[20]*pB[17] + pA[21]*pB[23] + pA[22]*pB[29] + pA[23]*pB[35];
  pOut[24] = pA[24]*pB[0] + pA[25]*pB[6] + pA[26]*pB[12] + pA[27]*pB[18] + pA[28]*pB[24] + pA[29]*pB[30];
  pOut[25] = pA[24]*pB[1] + pA[25]*pB[7] + pA[26]*pB[13] + pA[27]*pB[19] + pA[28]*pB[25] + pA[29]*pB[31];
  pOut[26] = pA[24]*pB[2] + pA[25]*pB[8] + pA[26]*pB[14] + pA[27]*pB[20] + pA[28]*pB[26] + pA[29]*pB[32];
  pOut[27] = pA[24]*pB[3] + pA[25]*pB[9] + pA[26]*pB[15] + pA[27]*pB[21] + pA[28]*pB[27] + pA[29]*pB[33];
  pOut[28] = pA[24]*pB[4] + pA[25]*pB[10] + pA[26]*pB[16] + pA[27]*pB[22] + pA[28]*pB[28] + pA[29]*pB[34];
  pOut[29] = pA[24]*pB[5] + pA[25]*pB[11] + pA[26]*pB[17] + pA[27]*pB[23] + pA[28]*pB[29] + pA[29]*pB[35];
  pOut[30] = pA[30]*pB[0] + pA[31]*pB[6] + pA[32]*pB[12] + pA[33]*pB[18] + pA[34]*pB[24] + pA[35]*pB[30];
  pOut[31] = pA[30]*pB[1] + pA[31]*pB[7] + pA[32]*pB[13] + pA[33]*pB[19] + pA[34]*pB[25] + pA[35]*pB[31];
  pOut[32] = pA[30]*pB[2] + pA[31]*pB[8] + pA[32]*pB[14] + pA[33]*pB[20] + pA[34]*pB[26] + pA[35]*pB[32];
  pOut[33] = pA[30]*pB[3] + pA[31]*pB[9] + pA[32]*pB[15] + pA[33]*pB[21] + pA[34]*pB[27] + pA[35]*pB[33];
  pOut[34] = pA[30]*pB[4] + pA[31]*pB[10] + pA[32]*pB[16] + pA[33]*pB[22] + pA[34]*pB[28] + pA[35]*pB[34];
  pOut[35] = pA[30]*pB[5] + pA[31]*pB[11] + pA[32]*pB[17] + pA[33]*pB[23] + pA[34]*pB[29] + pA[35]*pB[35];

  return out;
}

//////////////////////////////////////////////////
// This function is AUTOGENERATED. Do not edit! //
//////////////////////////////////////////////////
inline ActsMatrix<6, 6> transpose(const ActsMatrix<6, 6>& A) {
  ActsMatrix<6, 6> out;
  double* pOut = out.data();

  const double* pA = A.data();

  pOut[0] = pA[0];
  pOut[1] = pA[6];
  pOut[2] = pA[12];
  pOut[3] = pA[18];
  pOut[4] = pA[24];
  pOut[5] = pA[30];
  pOut[6] = pA[1];
  pOut[7] = pA[7];
  pOut[8] = pA[13];
  pOut[9] = pA[19];
  pOut[10] = pA[25];
  pOut[11] = pA[31];
  pOut[12] = pA[2];
  pOut[13] = pA[8];
  pOut[14] = pA[14];
  pOut[15] = pA[20];
  pOut[16] = pA[26];
  pOut[17] = pA[32];
  pOut[18] = pA[3];
  pOut[19] = pA[9];
  pOut[20] = pA[15];
  pOut[21] = pA[21];
  pOut[22] = pA[27];
  pOut[23] = pA[33];
  pOut[24] = pA[4];
  pOut[25] = pA[10];
  pOut[26] = pA[16];
  pOut[27] = pA[22];
  pOut[28] = pA[28];
  pOut[29] = pA[34];
  pOut[30] = pA[5];
  pOut[31] = pA[11];
  pOut[32] = pA[17];
  pOut[33] = pA[23];
  pOut[34] = pA[29];
  pOut[35] = pA[35];

  return out;
}
