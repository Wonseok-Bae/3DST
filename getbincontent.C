void getbincontent()
{
    TFile * file = new TFile("/home/local1/DUNE/hist_2D_purity.root"); //root파일 불러오기. 빨간색이 데이터

    TH1F * cut_arm = new TH1F; //1차원 그래프 변수 만들기
    cut_arm = (TH1F*)file->Get("cut_arm"); //1차원 그래프 변수에, 해당하는 그래프 데이터 집어넣기(변수가 그래프의 데이터를 포인터로 가리기케 하기) 빨간색이 데이터
    
    for(int i = 1; i < 27; i++) //이 루프동안
    {
        cout<<i<<","<<cut_arm->GetBinContent(i)<<endl;;//그래프 안의 데이터 꺼내기. 아래 왼쪽부터 1,2,3
    }

      
    TH2F * cut_2D = new TH2F; //2차원 그래프 변수 만들기
    cut_2D = (TH2F*)file->Get("cut_2D"); //2차원 그래프 변수에, 해당하는 그래프 데이터 집어넣기(변수가 그래프의 데이터를 포인터로 가리기케 하기) 빨간색이 데이터

    for(int i = 1; i < 100; i++) //이 루프동안
    {
        cout<<i<<","<<cut_2D->GetBinContent(i)<<endl;;//그래프 안의 데이터 꺼내기. 아래 왼쪽부터 1,2,3
    }
}
