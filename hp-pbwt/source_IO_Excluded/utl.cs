using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO.Compression;
namespace HP_PBWT
{
    public class utl
    {
        public static List<BitArray> RandomPanelGenerator(int nSite, int nHap)
        {
            Console.Write(DateTime.Now + " Generating Random Panel " + nSite + " by " + nHap + " ...");

            List<BitArray> panel = new List<BitArray>();
            for (int i = 0; i < nSite; i++)
            {
                panel.Add(new BitArray(nHap));

            }
            //Console.WriteLine(DateTime.Now + " Holder Created.");

            Parallel.For(0, nSite, (s) =>
            {
                Random rnd = new Random((int)DateTime.Now.Ticks + Thread.CurrentThread.ManagedThreadId);

                for (int h = 0; h < nHap; h++)
                {
                    if (rnd.Next(2) == 1)
                    {
                        panel[s][h] = true;
                    }
                    else
                    {
                        panel[s][h] = false;
                    }

                }

            });
            Console.WriteLine(" Done.");

            return panel;
        }
        

public static List<BitArray> LoadVCF(string filePath)
{
    Console.WriteLine(DateTime.Now + " Loading panel from VCF: " + filePath);

    StreamReader reader;
    if (filePath.EndsWith(".vcf.gz"))
    {
        var fileStream = File.OpenRead(filePath);
        var gzipStream = new GZipStream(fileStream, CompressionMode.Decompress);
        reader = new StreamReader(gzipStream);
    }
    else
    {
        reader = new StreamReader(filePath);
    }

    List<BitArray> panel = new List<BitArray>();
    int numSamples = 0;
    bool headerParsed = false;

    while (!reader.EndOfStream)
    {
        string line = reader.ReadLine();
        if (line.StartsWith("##")) continue;

        if (line.StartsWith("#CHROM"))
        {
            var tokens = line.Split('\t');
            numSamples = (tokens.Length - 9) * 2; // each sample has 2 haplotypes
            headerParsed = true;
            continue;
        }

        if (!headerParsed) continue;

        var fields = line.Split('\t');
        var genotypes = fields.Skip(9).ToArray();
        BitArray site = new BitArray(numSamples);

        for (int i = 0; i < genotypes.Length; i++)
        {
            if (genotypes[i].Length >= 3)
            {
                char a = genotypes[i][0];
                char b = genotypes[i][2];
                site[2 * i] = (a == '1');
                site[2 * i + 1] = (b == '1');
            }
        }

        panel.Add(site);
    }

    reader.Close();
    Console.WriteLine("VCF loaded with " + panel.Count + " sites and " + numSamples + " haplotypes.");
    return panel;
}

    }
}
