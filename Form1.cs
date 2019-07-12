using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Runtime.Serialization.Formatters.Binary;

namespace WindowsFormsApp3
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        Network test;
        string[] Right1, Wrong;

        private void Form1_Load(object sender, EventArgs e)
        {
            test = new Network(250,250,25,25,5,2,2);
        }

        private double[][,] ImageConvert(string File)
        {
            double[][,] Out = new double[3][,];
            Color col;
            Bitmap img =new Bitmap(Image.FromFile(openFileDialog1.FileName));
            Out[0] = new double[img.Width, img.Height];
            Out[1] = new double[img.Width, img.Height];
            Out[2] = new double[img.Width, img.Height];
            for (int i = 0; i < img.Height; i++)
                for (int j = 0; j < img.Width; j++)
                {
                    col = img.GetPixel(j, i);
                    Out[0][j, i] = col.R;
                    Out[1][j, i] = col.G;
                    Out[2][j, i] = col.B;
                }

            return Out;
        }

        private void button1_Click(object sender, EventArgs e)
        {
            this.BackColor = Color.Red;
            this.Update();

            openFileDialog1.ShowDialog();
            double[] Out = test.Compute(ImageConvert(openFileDialog1.FileName));

            string[] lines = new string[Out.Length];
            for (int i=0; i < Out.Length; i++)
            {
                lines[i] = Out[i].ToString();
            }
            richTextBox1.Lines = lines;

            this.BackColor = Color.Green;
            this.Update();
        }



        private void SaveNetwork(string file)
        {
            BinaryFormatter formatter = new BinaryFormatter();
            using (System.IO.FileStream stream = System.IO.File.OpenWrite(file))
            {
                formatter.Serialize(stream, test);
            }
        }

        private void LoadNetwork(string file)
        {
            BinaryFormatter formatter = new BinaryFormatter();
            using (System.IO.FileStream stream = System.IO.File.OpenRead(file))
            {
                test = (Network)formatter.Deserialize(stream);
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            this.BackColor = Color.Red;
            this.Update();

            var result = saveFileDialog1.ShowDialog();
            if (result==DialogResult.OK || result == DialogResult.Yes)
            SaveNetwork(saveFileDialog1.FileName);

            this.BackColor = Color.Green;
            this.Update();
        }

        private void button3_Click(object sender, EventArgs e)
        {
            this.BackColor = Color.Red;
            this.Update();

            var result = openFileDialog1.ShowDialog();
            if (result == DialogResult.OK || result == DialogResult.Yes)
                LoadNetwork(openFileDialog1.FileName);

            this.BackColor = Color.Green;
            this.Update();
        }

        private void button5_Click(object sender, EventArgs e)
        {
            openFileDialog1.Multiselect = true;
            openFileDialog1.ShowDialog();
            Right1 = openFileDialog1.FileNames;
            button5.BackColor = Color.Green;

            openFileDialog1.Multiselect = false;
        }

        private void button6_Click(object sender, EventArgs e)
        {
            openFileDialog1.Multiselect = true;
            openFileDialog1.ShowDialog();
            Wrong = openFileDialog1.FileNames;
            button6.BackColor = Color.Green;

            openFileDialog1.Multiselect = false;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            this.BackColor = Color.Red;
            this.Update();

            int x = Right1.Length;
            if (x > Wrong.Length)
                x = Wrong.Length;
            double[] a, b;
            a = new double[2] { 1,-1 };
            b = new double[2] { -1, 1 };
            for (int i = 0; i < x; i++)
            {
                test.Teach(ImageConvert(Right1[i]), a);
                test.Teach(ImageConvert(Wrong[i]), b);
            }

            this.BackColor = Color.Green;
            this.Update();
        }


    }
}
