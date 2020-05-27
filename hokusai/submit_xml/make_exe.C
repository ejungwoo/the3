void make_exe()
{
  auto numSplits = 20;

  ofstream file132("exe_132.sh");
  for (auto iSplit=0; iSplit<numSplits; ++iSplit)
    file132 << Form("root -q -b -l run_analysis_xml.C\\(132,%d,%d\\) > log/log_132_%d 2>&1",numSplits,iSplit,iSplit) << endl;

  ofstream file108("exe_108.sh");
  for (auto iSplit=0; iSplit<numSplits; ++iSplit)
    file108 << Form("root -q -b -l run_analysis_xml.C\\(108,%d,%d\\) > log/log_108_%d 2>&1 ",numSplits,iSplit,iSplit) << endl;
}
