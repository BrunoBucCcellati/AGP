#include "ALGORITHM.h"
namespace AGP {
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
		}
	protected:
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Button^ button1;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^ chart1;
	protected:
	private:
		System::ComponentModel::Container ^components;
#pragma region Windows Form Designer generated code
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataVisualization::Charting::ChartArea^ chartArea1 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^ legend1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^ series1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^ series2 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^ series3 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^ series4 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Title^ title1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Title());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->SuspendLayout();
			this->button1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button1->Location = System::Drawing::Point(1729, 12);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(201, 32);
			this->button1->TabIndex = 0;
			this->button1->Text = L"TEST";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			chartArea1->AxisX->Title = L"NUM OF ITERATIONS";
			chartArea1->AxisX->TitleFont = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			chartArea1->AxisY->Title = L"% CORRECT";
			chartArea1->AxisY->TitleFont = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			chartArea1->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea1);
			legend1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			legend1->IsTextAutoFit = false;
			legend1->Name = L"Legend1";
			this->chart1->Legends->Add(legend1);
			this->chart1->Location = System::Drawing::Point(12, 50);
			this->chart1->Name = L"chart1";
			series1->BorderWidth = 2;
			series1->ChartArea = L"ChartArea1";
			series1->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series1->Color = System::Drawing::Color::Blue;
			series1->Legend = L"Legend1";
			series1->LegendText = L"AGP_R2";
			series1->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Circle;
			series1->Name = L"Series1";
			series2->BorderWidth = 2;
			series2->ChartArea = L"ChartArea1";
			series2->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series2->Color = System::Drawing::Color::Fuchsia;
			series2->Legend = L"Legend1";
			series2->LegendText = L"AGP_LNA_R2";
			series2->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Cross;
			series2->Name = L"Series2";
			series3->BorderWidth = 2;
			series3->ChartArea = L"ChartArea1";
			series3->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series3->Legend = L"Legend1";
			series3->LegendText = L"AGP_R1";
			series3->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Triangle;
			series3->Name = L"Series3";
			series4->BorderWidth = 2;
			series4->ChartArea = L"ChartArea1";
			series4->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Line;
			series4->Legend = L"Legend1";
			series4->LegendText = L"AGP_LNA_R1";
			series4->MarkerStyle = System::Windows::Forms::DataVisualization::Charting::MarkerStyle::Diamond;
			series4->Name = L"Series4";
			this->chart1->Series->Add(series1);
			this->chart1->Series->Add(series2);
			this->chart1->Series->Add(series3);
			this->chart1->Series->Add(series4);
			this->chart1->Size = System::Drawing::Size(1918, 661);
			this->chart1->TabIndex = 1;
			this->chart1->Text = L"chart1";
			title1->Font = (gcnew System::Drawing::Font(L"Tahoma", 14.25F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			title1->Name = L"Title1";
			title1->Text = L"OPERATING CHARACTERISTICS OF AGP";
			this->chart1->Titles->Add(title1); 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1942, 723);
			this->Controls->Add(this->chart1);
			this->Controls->Add(this->button1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->ResumeLayout(false);
		}
#pragma endregion
	private:
		PeanoCurve_2D^ Curve = gcnew PeanoCurve_2D(11, List::Top, 0, 1, 0, 1);
		PeanoCurve_2D^ Curve_Minus_PI_Na_Dva = gcnew PeanoCurve_2D(11, List::Right, 0, 1, 0, 1);
		cli::array<unsigned short>^ procent_correct = gcnew cli::array<unsigned short>(1001);
		cli::array<unsigned short>^ procent_correct_LNA = gcnew cli::array<unsigned short>(1001);
		double correctness, MIN_GRSH = -9.79812, MIN_SH = -1.84034;
		unsigned short i = 1000;
		System::Void button1_Click(System::Object^ sender, System::EventArgs^ e)
		{
			chart1->Series[0]->Points->Clear();
			chart1->Series[1]->Points->Clear();
			chart1->Series[2]->Points->Clear();
			chart1->Series[3]->Points->Clear();
			while (i > 0)
			{
				procent_correct[i] = unsigned short();
				procent_correct_LNA[i--] = unsigned short();
			}
			chart1->ChartAreas[0]->AxisX->Minimum = 0;
			chart1->ChartAreas[0]->AxisX->Maximum = 1000;
			chart1->ChartAreas[0]->AxisY->Minimum = 0;
			chart1->ChartAreas[0]->AxisY->Maximum = 100;
			chart1->ChartAreas[0]->AxisX->MajorGrid->Interval = 50;
			chart1->ChartAreas[0]->AxisY->MajorGrid->Interval = 10;
			std::deque<double> Extr = Base_LNA_1_2_Mer_AGP(false, 2, 1, Curve, Curve_Minus_PI_Na_Dva, 1.1);
			std::deque<double> Extr_LNA = Base_LNA_1_2_Mer_AGP(true, 2, 1, Curve, Curve_Minus_PI_Na_Dva, 1.4);
			Extr.pop_front();
			Extr_LNA.pop_front();
			correctness = 0.001;
			while (Extr.empty() == false)
			{
				i++;
				if (Extr.front() - MIN_GRSH < correctness || Extr.front() < MIN_GRSH)
				{
					procent_correct[i] += 1 + procent_correct[i - 1];
				}
				else
				{
					procent_correct[i] += procent_correct[i - 1];
				}
				chart1->Series[0]->Points->AddXY(i, procent_correct[i] * 100.0 / i);
				Extr.pop_front();
				if (Extr_LNA.empty() == false)
				{
					if (Extr_LNA.front() - MIN_GRSH < correctness || Extr_LNA.front() < MIN_GRSH)
					{
						procent_correct_LNA[i] += 1 + procent_correct_LNA[i - 1];
					}
					else
					{
						procent_correct_LNA[i] += procent_correct_LNA[i - 1];
					}
					chart1->Series[1]->Points->AddXY(i, procent_correct_LNA[i] * 100.0 / i);
					Extr_LNA.pop_front();
				}
			}
			Extr = Base_LNA_1_2_Mer_AGP();
			Extr_LNA = Base_LNA_1_2_Mer_AGP(true);
			Extr.pop_front();
			Extr_LNA.pop_front();
			while (i > 0)
			{
				procent_correct[i] = unsigned short();
				procent_correct_LNA[i--] = unsigned short();
			}
			correctness = 0.0001;
			while (Extr.empty() == false)
			{
				i++;
				if (Extr.front() - MIN_SH < correctness || Extr.front() < MIN_SH)
				{
					procent_correct[i] += 1 + procent_correct[i - 1];
				}
				else
				{
					procent_correct[i] += procent_correct[i - 1];
				}
				chart1->Series[2]->Points->AddXY(i, procent_correct[i] * 100.0 / i);
				Extr.pop_front();
				if (Extr_LNA.empty() == false)
				{
					if (Extr_LNA.front() - MIN_SH < correctness || Extr_LNA.front() < MIN_SH)
					{
						procent_correct_LNA[i] += 1 + procent_correct_LNA[i - 1];
					}
					else
					{
						procent_correct_LNA[i] += procent_correct_LNA[i - 1];
					}
					chart1->Series[3]->Points->AddXY(i, procent_correct_LNA[i] * 100.0 / i);
					Extr_LNA.pop_front();
				}
			}
		}
	};
}