/*

Adjacency count

Given a sequence of 0 and 1's 
The adjacency count of the sequence is the maximum length of an interval
of equal values

Additional condition :
Exactly one of the values must be inverted

Fast algorithm :
Scan over values and keep track of adjacency sequence.
Notice that a seqence can be extended by allowing for exactly one value to be inverted 
if the value is of the opposite value of the members of the adjacency sequence
If during the scan an adjacency sequence can only be extended by inverting the current value for a second
time : the adjancency sequence ends.

*/


#include <vector>
#include <list>
#include <algorithm>
#include <iostream>


typedef std::vector<int> tIntVec;


// Test cases
tIntVec test_1 = { 1,1,1,0,1,0,0,0,1,1,1,0,0,1,1,0 };
tIntVec test_2 = { 1,1,1,0,1,0,0,0,1,1,1,1,1,0,0,1,1,0 };


// Find maximum adjacency count of a sequence
int
FindMaximumAdjacencyCount(const tIntVec& inSeq)
{
    if (inSeq.empty())
    {
        return 0;
    }

    int max_length(1);

    int length(1);
    int curr_val(inSeq[0]);
    for (auto s_it (inSeq.begin() + 1); s_it != inSeq.end(); ++s_it)
    {
        if (*s_it == curr_val)
        {
            ++length;            
        }
        else
        {
            max_length = std::max(max_length, length);
            curr_val = *s_it;
            length = 1;
        }
    }   

    // Compare with current adjacency sequence 
    return std::max(max_length, length);
}


// Find maximum adjacency count of a sequence where one value must be inverted
// Algorithm : iterate over complete sequence, invert value at iterator and run FindMaximumAdjacencyCount() over created sequence
int
FindMaximumAdjacencyCountWithOneElementInvertedSimple(const tIntVec& inSeq)
{
    tIntVec test_seq(inSeq);

    int max_seq(0);
    for (auto s_it(test_seq.begin()); s_it != test_seq.end(); ++s_it)
    {
        // Invert one item
        int val(*s_it);
        *s_it = (val == 0) ? 1 : 0;
        max_seq = std::max(max_seq, FindMaximumAdjacencyCount(test_seq));
        *s_it = val;
    }

    return max_seq;
}


// Adjacency sequence counter utility class
// 
class AdjSeqCounters
{
public:     
                                // Constructor
                                AdjSeqCounters(int inTrackValue);

                                // Update counters
                                // Return adjacency count if sequence ends
    void                        Update(int inValue);

                                // Return adjacency count for valid sequences 
                                // Takes into account that one value must be inverted
    int                         GetAdjacencySequenceLengthAtEnd();

private:
                                // Value of individual elements for which the counters are tracking
    int                         mTrackValue;

                                // Current length of adjacency sequence for which not yet one value is inverted
    int                         mSeqCount;

                                // Current length of adjacency sequence for which exactly one value is inverted
    int                         mSeqCountOneInv;

                                // Maximum length of adjacency sequence with exactly one value inverted
    int                         mMaxSeqCount;
};


// Constructor
AdjSeqCounters::AdjSeqCounters(int inTrackValue) :
    mTrackValue     (inTrackValue),
    mSeqCount       (0),
    mSeqCountOneInv (0),
    mMaxSeqCount    (0)
{
}

// Update counters
// Track maximum adjacency count if adjacency sequence ends
void
AdjSeqCounters::Update(int inValue)
{
    if (inValue == mTrackValue)
    {
        // Current value is equal to value of tracked sequences
        ++mSeqCount;

        if (mSeqCountOneInv > 0)
        {
            // Sequence with one inverted value exists
            ++mSeqCountOneInv;
        }
    }
    else
    {
        // Current value is opposite of value of tracked sequences
        // Sequence with one inverted value ends
        mMaxSeqCount = std::max(mMaxSeqCount, mSeqCountOneInv);
        mSeqCountOneInv = 0;
        if (mSeqCount > 0)
        {
            // Sequence with no inverted values is extended with one inverted value
            mSeqCountOneInv = mSeqCount + 1;

            // Sequence with no inverted value is reset
            mSeqCount = 0;
        }
    }
}


// Return adjacency count for valid sequences 
// Takes into account that one value must be inverted
int             
AdjSeqCounters::GetAdjacencySequenceLengthAtEnd()
{
    return std::max(
        std::max(
            mSeqCount - 1, 
            mSeqCountOneInv),
        mMaxSeqCount);
}


// Find maximum adjacency count of a sequence with one value inverted
// In one scan over the sequence
int
FindMaximumAdjacencyCountWithOneElementInvertedFast(const tIntVec& inSeq)
{
    AdjSeqCounters seq_counters_zeros(0);
    AdjSeqCounters seq_counters_ones(1);

    int max_seq(0);
    for (auto s_it(inSeq.begin()); s_it != inSeq.end(); ++s_it)
    {
        seq_counters_zeros.Update(*s_it);
        seq_counters_ones.Update(*s_it);
    }
    
    return std::max(
        std::max(
            seq_counters_zeros.GetAdjacencySequenceLengthAtEnd(),
            seq_counters_ones.GetAdjacencySequenceLengthAtEnd()
        ),
        max_seq);
}


// Test parameters

// Number of tests
const int TEST_COUNT     =  100;

// Maximum length of a sequence
const int MAX_SEQ_LENGTH = 1000;


int
main()
{
    // Human verifiable test cases
    std::cout << "Max adjacency count test_1 : " << FindMaximumAdjacencyCount(test_1) << std::endl;
    std::cout << "Max adjacency count test_2 : " << FindMaximumAdjacencyCount(test_2) << std::endl;

    std::cout << "Max one flipped seq test_1 simple : " << FindMaximumAdjacencyCountWithOneElementInvertedSimple(test_1) << std::endl;
    std::cout << "Max one flipped seq test_2 simple : " << FindMaximumAdjacencyCountWithOneElementInvertedSimple(test_2) << std::endl;

    std::cout << "Max one flipped seq test_1 fast : " << FindMaximumAdjacencyCountWithOneElementInvertedFast(test_1) << std::endl;
    std::cout << "Max one flipped seq test_2 fast : " << FindMaximumAdjacencyCountWithOneElementInvertedFast(test_2) << std::endl;

    // Computer generated test cases
    // Verified by comparison with simple algorithm
    int test_case_fail_count(0);
    for (int t_ix(0); t_ix < TEST_COUNT; ++t_ix)
    {
        // Generate testcase
        tIntVec seq(rand() % MAX_SEQ_LENGTH);
        for_each(seq.begin(), seq.end(), [](int &outVal) { outVal = rand() % 2; });

        // Test
        int max_adj_count_simpe(FindMaximumAdjacencyCountWithOneElementInvertedSimple(seq));
        int max_adj_count_fast(FindMaximumAdjacencyCountWithOneElementInvertedFast(seq));

        if (max_adj_count_simpe != max_adj_count_fast)
        {
            ++test_case_fail_count;
        }
    }

    std::cout << "Count of test cases for which a difference was found between simple and fast algorithm : " << test_case_fail_count << std::endl;

    return 0;
}
