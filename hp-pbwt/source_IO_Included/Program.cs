using HP_PBWT;
using System.Diagnostics;

public class Program
{
    #region debug use
    public static StreamReader dbgSR;

    #endregion

    public static List<int> LM_Collection = new List<int>();
    public static int nHap = 0;
    public static int nIndv = 0;
    public static List<char[]> BufferArr = new List<char[]>();
    public static bool semiRun = false;

    public static PBWT.Pal pal;


    public static utl.BufferWriter BW;

    public static int THD_nHapSemiRun = 10000;
    public static int rtrBlockSize = 500000000;
    //500000
    public static int rtrN_Block = 20;

    #region user inputs
    public static int nThread = 0;
    public static int LLM_Len = 0;

    public static int LC_THD = 0;
    #endregion


    private static void Main(string[] args)
    {


        if (args.Count() != 4 && args.Count() != 7)
        {
            Console.WriteLine("HP-PBWT");
            Console.WriteLine("<Input VCF path,string> <Output path,string> <# of Thread,int> <L-Long Match Length,int. <1 for SetMaxMatch.>");
            Console.WriteLine("<Input VCF path,string> <Output path,string> <# of Thread,int> <L-Long Match Length,int. <1 for SetMaxMatch.> --SetBuffer <Buffer Size> <# of Buffer>"); ;
            return;
        }

        string inVCF, outFile;
        inVCF = args[0];
        outFile = args[1];
        nThread = Convert.ToInt32(args[2]);
        LLM_Len = Convert.ToInt32(args[3]);

        Console.WriteLine("HP-PBWT");
        Console.WriteLine("VCF: " + inVCF);
        Console.WriteLine("Output: " + outFile);
        Console.WriteLine("nThread: " + nThread);
        Console.WriteLine("Length: " + LLM_Len);

        if (args.Count() == 7)
        {
            rtrBlockSize = Convert.ToInt32(args[5]);
            rtrN_Block = Convert.ToInt32(args[6]);
        }



        LC_THD = LLM_Len - 1;


        Stopwatch spw = Stopwatch.StartNew();
        spw.Start();

        Run(inVCF, outFile);

        spw.Stop();
        Console.WriteLine();
        Console.WriteLine("MS: " + spw.ElapsedMilliseconds);

    }


    public static void Run(string VCF_Path, string outPath)
    {
        utl.RoundTableReaderV2 rtr = new utl.RoundTableReaderV2(VCF_Path);
        BufferArr = rtr.BufferArr;

        BW = new utl.BufferWriter(outPath);


        Task t1 = Task.Run(() => rtr.BlockReader());

        rtr.ProbTop();


        Program.pal = new PBWT.Pal();


        Task t2 = Task.Run(() => rtr.LineParser());

        Task t3;
        if (LLM_Len < 1)
        {
            t3 = Task.Run(() => rtr.LineTaker_SMM());
        }
        else
        {
            t3 = Task.Run(() => rtr.LineTaker_LLM());
        }

        Task t4 = Task.Run(() => BW.Run());

        Task.WaitAll(t1, t2, t3, t4);


    }
}