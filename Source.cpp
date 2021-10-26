/*
Pligiarism Detection using string matching algorithms (Naive - KMP - Boyer Moore - Rabin Karp)
*/

#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<unordered_map>
#include<Windows.h>
using namespace std;

# define d 256 //number of possible chars, will be used in both Boyer Moore and Rabin Karp

//Approximate String Matching (BruteForce) algorithm
int ApproximateStringMatching_Search(string pat, string txt)
{
	int M = pat.size();
	int N = txt.size();

	for (int i = 0; i <= N - M; i++) {
		int distance = 0;
		for (int j = 0; j < M; j++)
			if (pat[j] != txt[i + j])
				distance++;
		if (distance == 0)
			return i;
	}
	return -1;
}

//KMP algorithm, the two coming functions
// Calculating the longest proper suffix/prefix array for the given pattern 
void calcLPS(string pat, int L, int* lps)
{
	// indicates the length of the previous LPS
	int ln = 0;

	lps[0] = 0; // Always 0 because there is no lps for a single character

	//computing lps[i]
	int i = 1;
	while (i < L)
	{
		if (pat[i] == pat[ln])
		{
			ln++;
			lps[i] = ln;
			i++;
		}

		else //not equal
		{
			if (ln == 0)
			{
				lps[i] = 0;
				i++;
			}
			else
				ln = lps[ln - 1];
		}
	}
}

int KMP_Algorithm(string pat, string txt)
{
	int L = pat.length();
	int K = txt.length();

	//array to store the longest proper/perfix suffix of the given pattern
	int* lps = new int[L];

	//calling the function to compute the lps array
	calcLPS(pat, L, lps);

	int i = 0;//txt 
	int j = 0;//pat

	while (i < K)
	{
		if (pat[j] == txt[i]) // if equal, means the two chars are identical, so continue rolling over pat and txt
		{
			j++;
			i++;
		}

		if (j == L) //meaning that we have rolled over the whole pattern, so it's found inside the txt
			return i - j;   //the start index of matching. its value (>=0) means that the pattern is found

		// a mismatch is detected
		else if (i < K && pat[j] != txt[i])
		{

			if (j == 0)
				i = i + 1;//only role over txt as the current index over pat will be at its very beginning
			else
				j = lps[j - 1];//using the lps array to avoid rematching characters that will certainly match
		}
	}
	return -1; //meaning that the pattern is not found
}

//Boyer Moore algorithm, the two coming functions
void badCharHeuristic(string str, int size, int badchar[d])
{
	int i;

	// Initialize all occurrences as -1
	for (i = 0; i < d; i++)
		badchar[i] = -1;

	// Fill the actual value of last occurrence
	// of a character
	for (i = 0; i < size; i++)
		badchar[(int)str[i]] = i;
}

/* A pattern searching function that uses Bad
Character Heuristic of Boyer Moore Algorithm */
int Boyersearch(string pat, string txt)
{
	int m = pat.size();
	int n = txt.size();

	int badchar[d];

	/* Fill the bad character array by calling
	the preprocessing function badCharHeuristic()
	for given pattern */
	badCharHeuristic(pat, m, badchar);

	int s = 0; // s is shift of the pattern with
				// respect to text
	while (s <= (n - m))
	{
		int j = m - 1;

		/* Keep reducing index j of pattern while
		characters of pattern and text are
		matching at this shift s */
		while (j >= 0 && pat[j] == txt[s + j])
			j--;

		/* If the pattern is present at current
		shift, then index j will become -1 after
		the above loop */
		if (j < 0)
		{
			return s;
			/* Shift the pattern so that the next
			character in text aligns with the last
			occurrence of it in pattern.
			The condition s+m < n is necessary for
			the case when pattern occurs at the end
			of text */
			s += (s + m < n) ? m - badchar[txt[s + m]] : 1;
		}

		else
			/* Shift the pattern so that the bad character
			in text aligns with the last occurrence of
			it in pattern. The max function is used to
			make sure that we get a positive shift.
			We may get a negative shift if the last
			occurrence of bad character in pattern
			is on the right side of the current
			character. */
			s += max(1, j - badchar[txt[s + j]]);
	}
	return -1;
}

//Rabin Karp Algorithm
int Rabin_Karp_Algorithm(string pat, string txt, int pri)
{
	int L = pat.length();
	int K = txt.length();
	int i, j;
	int pat_hash = 0; // hash value of pat
	int txt_hash = 0; // hash value of txt
	int p = 1;

	// final value of p = pow(d, L-1)%pri
	for (i = 0; i < L - 1; i++)
		p = (p * d) % pri;

	// computing the hash value of pat and first window of t
	for (i = 0; i < L; i++)
	{
		pat_hash = (d * pat_hash + pat[i]) % pri;
		txt_hash = (d * txt_hash + txt[i]) % pri;
	}

	// role over txt by the size of pat one by one
	for (i = 0; i <= K - L; i++)
	{
		if (pat_hash == txt_hash)//if hash values are equal, then the characters're potentially identical
		{
			//check character by character
			for (j = 0; j < L; j++)
			{
				if (txt[i + j] != pat[j])//if a mistach is found, break as that indicates that they are not identical
					break;
			}

			if (j == L) //meaning that we have rolled over the whole pattern, so it's found inside the txt
				return i;  //the start index of matching. its value (>=0) means that the pattern is found
		}

		/* computing the next window hash value like we are rolling over txt; delete the value of the far left char and add
		the value of the far right char */
		if (i < K - L)
		{
			//using %pri is to make sure that the value of txt_hash will not exceed the max value to be stored in the variable
			txt_hash = (d * (txt_hash - txt[i] * p) + txt[i + L]) % pri;

			//converting -ve value to +ve if found
			if (txt_hash < 0)
				txt_hash = (txt_hash + pri);
		}
	}
	return -1; //meaning that the pattern is not found
}

int main()
{
	int option1, option2, option3;
	int index;
	int pri = 101; //prime number to be used in the Rabin Karp algoritm

	ifstream input;

	string t;
	string  pat, line;

	vector <string> p;	       // store the sentences of the potentially plagiarized doc seperately
	vector<string> textnames;  // store sources' names
	vector<string> texts;      // store the text inside the sources

	unordered_map<string, vector<int>> pligiarized_from; //store the text name of the document plagiarized from
	unordered_map<string, vector<int>> indexes;          //store the index at which plagiarism starts from in the document plagiarized from


	cout << "what do you want to do?\n1. check for pligiarism\n2. Exit\nEnter your choice:  ";
	cin >> option1;

	if (option1 == 1) {

		cout << "Type in the number of original texts to compare with: ";
		cin >> option3;

		cout << "Type in their names with their extensions: ";
		for (int i = 0; i < option3; i++)
		{
			cin >> t;
			textnames.push_back(t);
		}

		cout << "Type in the name of the potentially plagiarized document: ";
		cin >> t;
		textnames.push_back(t);

		//reading the text inside each source textfile 
		for (int i = 0; i < option3; i++)
		{
			input.open(textnames[i]);
			while (getline(input, line))
			{
				texts.push_back("");
				texts[i] += line + " ";
			}
			input.close();
		}

		//reading sentences of the potentially plagiarized doc
		input.open(textnames.back());
		while (getline(input, pat, '.'))
		{
			if (pat.front() == ' ')
				pat.erase(0, 1);
			p.push_back(pat);
		}
		input.close();

		//choosing an algorithm to be excuted
		cout << "Type in the number of the desired algorithm:\n1. Naive\n2. KMP\n3. Rabin Karp\n4. Boyer Moore\n5. Exit\nNumber: ";
		cin >> option2;

		if (option2 > 0 && option2 < 5)
		{
			//Approximate String Matching 
			if (option2 == 1)
			{
				for (int i = 0; i < p.size(); i++) //looping over each sentence of the test doc
				{
					for (int j = 0; j < option3; j++)  //looping over every source textfile
					{
						index = ApproximateStringMatching_Search(p[i], texts[j]);   //calling the function
						if (index != -1)                //if index!=-1 then it is plagiarized and index represents the starting index
						{
							pligiarized_from[textnames[j]].push_back(i + 1);  // adding the source textfile name
							indexes[textnames[j]].push_back(index);           // adding the start index of plagiarism
						}
					}
				}
			}


			//KMP
			else if (option2 == 2)
			{
				for (int i = 0; i < p.size(); i++) //looping over each sentence of the test doc
				{
					for (int j = 0; j < option3; j++) //looping over every source textfile
					{
						index = KMP_Algorithm(p[i], texts[j]);  //calling the function
						if (index != -1)              //if index!=-1 then it is plagiarized and index represents the starting index
						{
							pligiarized_from[textnames[j]].push_back(i + 1);  // adding the source textfile name
							indexes[textnames[j]].push_back(index);           // adding the start index of plagiarism
						}
					}
				}
			}

			//Rabin-Karp
			else if (option2 == 3)
			{
				for (int i = 0; i < p.size(); i++) //looping over each sentence of the test doc
				{
					for (int j = 0; j < option3; j++) //looping over every source textfile
					{
						index = Rabin_Karp_Algorithm(p[i], texts[j], pri);  //calling the function
						if (index != -1)          //if index!=-1 then it is plagiarized and index represents the starting index
						{
							pligiarized_from[textnames[j]].push_back(i + 1);  // adding the source textfile name
							indexes[textnames[j]].push_back(index);           // adding the start index of plagiarism
						}
					}
				}

			}
			//Boyer Moore
			else if (option2 == 4)
			{
				for (int i = 0; i < p.size(); i++) //looping over each sentence of the test doc
				{
					for (int j = 0; j < option3; j++) //looping over every source textfile
					{
						index = Boyersearch(p[i], texts[j]);  //calling the function
						if (index != -1)          //if index!=-1 then it is plagiarized and index represents the starting index
						{
							pligiarized_from[textnames[j]].push_back(i + 1);  // adding the source textfile name
							indexes[textnames[j]].push_back(index);           // adding the start index of plagiarism
						}
					}
				}

			}

			string name;

			if (pligiarized_from.size() == 0) //checking if plagiarized or not
				cout << "\nNot pligiarized" << endl;

			else //then, plagiarized
			{
				cout << "\nDocuments Pligiarized from:\n";

				for (auto& e : pligiarized_from)//loop over all added source textfiles to the plagiarized from list, then display them
				{
					name = e.first;  //know its name
					cout << name;

					//these are the ordinal number of the pligiarized statement in the checked document
					for (int i = 0; i < pligiarized_from[name].size(); i++)
					{
						cout << "\n- Statement " << pligiarized_from[name][i] << " in the test document ";
						cout << "and starting from index " << indexes[name][i] << " in the original text";
					}
					cout << "\n\n";
				}
			}
		}
	}
	return 0;
}