 
using System.Collections;
using System.Diagnostics;
using HP_PBWT;

public class Program
{
    #region user inputs


    public static int nSite = 0;
    public static int nHap = 0;
    public static int nThread = 0;
    public static int LLM_Len = 0;
    public static int nRun = 0;

    #endregion

    public static int MS_Cooldown = 60000;
    public static int LC_THD;
    public static List<BitArray> panel;
    public static int[] Sums;

    private static void Main(string[] args)
{
    if (args.Length != 5)
    {
        Console.WriteLine("Cách sử dụng:");
        Console.WriteLine("  dotnet run --project HP-PBWT.csproj -- <đường_dẫn_vcf> <số_luồng> <độ_dài_LLM> <số_lần_chạy> <thời_gian_cooldown_ms>");
        Console.WriteLine("Ví dụ:");
        Console.WriteLine("  dotnet run --project HP-PBWT.csproj -- data/example.vcf 4 4 1 0");
        return;
    }

    // 1) Đường dẫn file VCF
    string vcfPath = args[0];

    // 2) Số luồng
    nThread     = Convert.ToInt32(args[1]);

    // 3) Độ dài L-long match
    LLM_Len     = Convert.ToInt32(args[2]);

    // 4) Số lần chạy
    nRun        = Convert.ToInt32(args[3]);

    // 5) Thời gian ngủ giữa các lần benchmark (ms)
    MS_Cooldown = Convert.ToInt32(args[4]);

    // Thiết lập ngưỡng cho Pal
    LC_THD      = LLM_Len - 1;

    // Load panel từ VCF
    panel = utl.LoadVCF(vcfPath);

    // Cập nhật số site và số haplotypes
    nSite = panel.Count;
    nHap  = (panel.Count > 0 ? panel[0].Length : 0);

    Console.WriteLine("=== HP-PBWT Benchmark từ file VCF ===");
    Console.WriteLine($"Đường dẫn VCF   : {vcfPath}");
    Console.WriteLine($"Số site        : {nSite}");
    Console.WriteLine($"Số haplotypes  : {nHap}");
    Console.WriteLine($"Số luồng       : {nThread}");
    Console.WriteLine($"Độ dài LLM     : {LLM_Len}");
    Console.WriteLine($"Số lần chạy    : {nRun}");
    Console.WriteLine($"Cooldown       : {MS_Cooldown} ms");
    Console.WriteLine();

    // Mở log
    // using var swLog = new StreamWriter($"{DateTime.Now.Ticks}.log") { NewLine = "\n" };
     string logFile = Path.Combine(Directory.GetCurrentDirectory(), "hp_pbwt_result.log");
    using var swLog = new StreamWriter(logFile, append: false) { NewLine = "\n" };




    // Vòng lặp benchmark
    for (int r = 0; r < nRun; r++)
    {
        Sums = new int[nHap];
        Console.WriteLine($"{DateTime.Now}: Đang chờ cooldown {MS_Cooldown} ms...");
        Thread.Sleep(MS_Cooldown);

        Stopwatch spw = Stopwatch.StartNew();
        if (nThread != 1)
        {
            Console.WriteLine("Chạy HP-PBWT (song song)...");
            new PBWT.Pal().Run();
        }
        else
        {
            Console.WriteLine("Chạy PBWT tuần tự...");
            new PBWT.Seq().Run();
        }
        spw.Stop();

        // Ghi log và in kết quả
        swLog.WriteLine($"{nSite}\t{nHap}\t{nThread}\t{LLM_Len}\t{MS_Cooldown}\t{spw.ElapsedMilliseconds}");
        Console.WriteLine($"Thời gian thực thi: {spw.ElapsedMilliseconds} ms\n");
    }
}





    public static void checkSum()
    {
        bool good = false;
        Parallel.For(0, nHap, (i, stat) =>
        {
            if (Program.Sums[i] > 0)
            {
                good = true;
                Console.WriteLine("Good Number. " + i + "  " + Program.Sums[i]);
                stat.Break();
            }


        });

        if (!good)
        {
            Console.WriteLine("XXXXX   ----All 0s----   XXXXX");
        }
    }

}