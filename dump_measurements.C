void dump_measurements() {
  TTree *t = (TTree*)gDirectory->Get("measurements");
  if (!t) {
    TFile *f = TFile::Open("measurements.root");
    t = (TTree*)f->Get("measurements");
  }

  FILE *out = fopen("dump.csv", "w");
  if (!out) {
    printf("Cannot open dump.csv for writing\n");
    return;
  }

  t->Draw("volume_id:layer_id:surface_id:extra_id", "", "goff");
  int n = t->GetSelectedRows();
  std::vector<int> vol(n), lay(n), sen(n), ext(n);
  for (int i = 0; i < n; i++) {
    vol[i] = (int)t->GetV1()[i];
    lay[i] = (int)t->GetV2()[i];
    sen[i] = (int)t->GetV3()[i];
    ext[i] = (int)t->GetV4()[i];
  }

  t->Draw("true_x:true_y:true_z:true_loc0", "", "goff");
  for (int i = 0; i < n; i++) {
    fprintf(out, "%d,%d,%d,%d,%f,%f,%f,%f\n",
            vol[i], lay[i], sen[i], ext[i],
            t->GetV1()[i], t->GetV2()[i], t->GetV3()[i], t->GetV4()[i]);
  }

  fclose(out);
  printf("Wrote %d entries to dump.csv\n", n);
}

