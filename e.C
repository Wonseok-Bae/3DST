void e()
{
    TFile * fi = new TFile("/Users/gwon/eff2D/threshold_vs_eff.root");
        TH2F * histo1 = new TH2F;
        TH2F * histo2 = new TH2F;
        TH2F * histo3 = new TH2F;
        TH2F * histo4 = new TH2F;
        TH2F * histo5 = new TH2F;
        TH2F * histo6 = new TH2F;
        TH2F * histo7 = new TH2F;
        TH2F * histo8 = new TH2F;
        TH2F * histo9 = new TH2F;
        TH2F * histo10 = new TH2F;

        TH2F * hhisto1 = new TH2F("test1", "test1;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto2 = new TH2F("test2", "test2;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto3 = new TH2F("test3", "test3;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto4 = new TH2F("test4", "test4;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto5 = new TH2F("test5", "test5;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto6 = new TH2F("test6", "test6;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto7 = new TH2F("test7", "test7;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto8 = new TH2F("test8", "test8;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto9 = new TH2F("test9", "test9;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        TH2F * hhisto10 = new TH2F("test10", "test10;Lever Arm [cm]; Time [ns]", 20, 0, 200, 25, 0, 25);
        
    histo1 = (TH2F*)fi->Get("efficiency_2D_1");
    histo2 = (TH2F*)fi->Get("efficiency_2D_2");
    histo3 = (TH2F*)fi->Get("efficiency_2D_3");
    histo4 = (TH2F*)fi->Get("efficiency_2D_4");
    histo5 = (TH2F*)fi->Get("efficiency_2D_5");
    histo6 = (TH2F*)fi->Get("efficiency_2D_6");
    histo7 = (TH2F*)fi->Get("efficiency_2D_7");
    histo8 = (TH2F*)fi->Get("efficiency_2D_8");
    histo9 = (TH2F*)fi->Get("efficiency_2D_9");
    histo10 = (TH2F*)fi->Get("efficiency_2D_10");

    for(int i = 1; i < 1000; i++)
    {
        if(histo1->GetBinContent(i) > 0)
            hhisto1->SetBinContent(i,histo1->GetBinContent(i));
        if(histo2->GetBinContent(i) > 0)
            hhisto2->SetBinContent(i,histo2->GetBinContent(i));
        if(histo3->GetBinContent(i) > 0)
            hhisto3->SetBinContent(i,histo3->GetBinContent(i));
        if(histo4->GetBinContent(i) > 0)
            hhisto4->SetBinContent(i,histo4->GetBinContent(i));
        if(histo5->GetBinContent(i) > 0)
            hhisto5->SetBinContent(i,histo5->GetBinContent(i));
        if(histo6->GetBinContent(i) > 0)
            hhisto6->SetBinContent(i,histo6->GetBinContent(i));
        if(histo7->GetBinContent(i) > 0)
            hhisto7->SetBinContent(i,histo7->GetBinContent(i));
        if(histo8->GetBinContent(i) > 0)
            hhisto8->SetBinContent(i,histo8->GetBinContent(i));
        if(histo9->GetBinContent(i) > 0)
            hhisto9->SetBinContent(i,histo9->GetBinContent(i));
        if(histo10->GetBinContent(i) > 0)
            hhisto10->SetBinContent(i,histo10->GetBinContent(i));
    }
    hhisto1->SetStats(false);
    hhisto2->SetStats(false);
    hhisto3->SetStats(false);
    hhisto4->SetStats(false);
    hhisto5->SetStats(false);
    hhisto6->SetStats(false);
    hhisto7->SetStats(false);
    hhisto8->SetStats(false);
    hhisto9->SetStats(false);
    hhisto10->SetStats(false);

    hhisto1->SetTitle("threshold : 0.1MeV");
    hhisto2->SetTitle("threshold : 0.2MeV");
    hhisto3->SetTitle("threshold : 0.3MeV");
    hhisto4->SetTitle("threshold : 0.4MeV");
    hhisto5->SetTitle("threshold : 0.5MeV");
    hhisto6->SetTitle("threshold : 0.6MeV");
    hhisto7->SetTitle("threshold : 0.7MeV");
    hhisto8->SetTitle("threshold : 0.8MeV");
    hhisto9->SetTitle("threshold : 0.9MeV");
    hhisto10->SetTitle("threshold : 1MeV");
    
    TCanvas * can = new TCanvas;
    can->Divide(3,4);
        can->cd(1);
        hhisto1->Draw("colz");
        can->cd(2);
        hhisto2->Draw("colz");
        can->cd(3);
        hhisto3->Draw("colz");
        can->cd(4);
        hhisto4->Draw("colz");
        can->cd(5);
        hhisto5->Draw("colz");
        can->cd(6);
        hhisto6->Draw("colz");
        can->cd(7);
        hhisto7->Draw("colz");
        can->cd(8);
        hhisto8->Draw("colz");
        can->cd(9);
        hhisto9->Draw("colz");
        can->cd(10);
        hhisto10->Draw("colz");
}
