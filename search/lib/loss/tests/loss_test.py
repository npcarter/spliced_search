import torch
from .. import loss

import numpy as np
def test_perfect_prediction():

    #create logits for extremely-good prediction
    true_1 = np.full(15, -100000.0)
    true_1 [0] = 100.0 #intergenic; large but finite logit so softmax/cross_entropy stays well-defined

    true_2 = np.full(15, -100000.0)
    true_2 [1] = 100.0 #start

    true_3 = np.full(15, -100000.0)
    true_3 [4] = 100.0 #exon1

    true_4 = np.full(15, -100000.0)
    true_4 [2] = 100.0 #stop

    true_5 = np.full(15, -100000.0)
    true_5 [0] = 100.0 #intergenic

    perfect_prediction  = torch.column_stack((torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])), 
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))))
    true_1 = np.zeros(15, dtype=float)
    true_1 [0] = 1.0 #intergenic

    true_2 = np.zeros(15, dtype=float)
    true_2 [1] = 1.0 #start

    true_3 = np.zeros(15, dtype=float)
    true_3 [4] = 1.0 #exon1

    true_4 = np.zeros(15, dtype=float)
    true_4 [2] = 1.0 #stop

    true_5 = np.zeros(15, dtype=float)
    true_5 [0] = 1.0 #intergenic

    ground_truth = torch.column_stack((torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])), 
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))))

    loss_model = loss.TiberiusLoss(use_f1=False, f1_factor=2.0)

    test_loss = loss_model.forward(perfect_prediction, ground_truth)
    assert (test_loss > -0.01 and test_loss < 0.01), f"Loss on perfect match was unexpected value {test_loss}"

    #retry with F1 loss turned on
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)
    
    test_loss = loss_model.forward(perfect_prediction, ground_truth)
        
    assert (test_loss > -0.01 and test_loss < 0.01), f"Loss on perfect match with F1 loss turned on was unexpected value {test_loss}"

    print("All tests passed for test_perfect_prediction()")


def test_awful_prediction():
    #create logits for extremely-strong prediction
    true_1 = np.full(15, -100000.0)
    true_1 [0] = 100000.0 #intergenic; large but finite logit so softmax/cross_entropy stays well-defined

    true_2 = np.full(15, -100000.0)
    true_2 [1] = 100000.0 #start

    true_3 = np.full(15, -100000.0)
    true_3 [4] = 100000.0 #exon1

    true_4 = np.full(15, -100000.0)
    true_4 [2] = 100000.0 #stop

    true_5 = np.full(15, -100000.0)
    true_5 [0] = 100000.0 #intergenic

    awful_prediction  = torch.column_stack((torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])), 
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))
                                       ))

    #create one-hot of completely-different ground truth
    true_1 = np.zeros(15, dtype=float)
    true_1 [3] = 1.0

    true_2 = np.zeros(15, dtype=float)
    true_2 [12] = 1.0 

    true_3 = np.zeros(15, dtype=float)
    true_3 [7] = 1.0

    true_4 = np.zeros(15, dtype=float)
    true_4 [14] = 1.0

    true_5 = np.zeros(15, dtype=float)
    true_5 [8] = 1.0 

    ground_truth = torch.column_stack((torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])), 
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5])),
                                       torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))
                                       ))

    loss_model = loss.TiberiusLoss(use_f1=False, f1_factor=2.0)

    test_loss = loss_model.forward(awful_prediction, ground_truth)

    #Here, we assume that Torch's CCE function works correctly
    assert (test_loss > 119999.0 and test_loss < 120000.0), f"Loss on awful prediction was unexpected value {test_loss}"


    #retry with f1_score turned on
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)

    test_loss1 = loss_model.forward(awful_prediction, ground_truth)

    #Here, we assume that Torch's CCE function works correctly
    assert (test_loss1 > 119999.0 and test_loss1 < 120000.0), f"Loss on awful prediction was unexpected value {test_loss}"
    assert test_loss1 > test_loss, "Adding F1 score didn't increase loss in awful prediction test"
    print("All tests passed for test_awful_prediction()")


def test_fp_penalty():
    #create one-hot of ground truth with no coding positions
    true_1 = np.zeros(15, dtype=float)
    true_1 [0] = 1.0 #intergenic
        
    true_2 = np.zeros(15, dtype=float)
    true_2 [6] = 1.0 #intron 
        
    true_3 = np.zeros(15, dtype=float)
    true_3 [7] = 1.0 #intron
        
    true_4 = np.zeros(15, dtype=float)
    true_4 [8] = 1.0 #intron
        
    true_5 = np.zeros(15, dtype=float)
    true_5 [0] = 1.0 #intergenic
        
    ground_truth = torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))
    #create logits for extremely-strong prediction
    true_1 = np.full(15, -100000.0)
    true_1 [0] = 100000.0 #intergenic; large but finite logit so softmax/cross_entropy stays well-defined

    true_2 = np.full(15, -100000.0)
    true_2 [1] = 100000.0 #start

    true_3 = np.full(15, -100000.0)
    true_3 [4] = 100000.0 #exon1

    true_4 = np.full(15, -100000.0)
    true_4 [2] = 100000.0 #stop

    true_5 = np.full(15, -100000.0)
    true_5 [0] = 100000.0 #intergenic

    fp_prediction  = torch.tensor(np.array([true_1, true_2, true_3, true_4, true_5]))

    #retry with f1_score turned on
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)

    test_loss1 = loss_model.forward(fp_prediction, ground_truth)

    #Here, we assume that Torch's CCE function works correctly
    assert (test_loss1 > 119999.9 and test_loss1 < 120000.6), f"Loss on FPR test was unexpected value {test_loss1}"
    print("All tests passed for test_fp_penalty()")

def test_state_reduce():
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)
    input_vector = torch.tensor([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 1, 0, 0, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 1, 0, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1,0,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1,0,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,1,0,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,1,0,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,1,0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,1]
                                   ])
    
    out = loss_model._reduce_states(input_vector)
    assert torch.equal(out, torch.tensor([0,2,4,2,3,4,1,1,1,2,3,4,2,3,4])), "Reduced state vector did not match expected"

def test_logit_reduce():
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)
    input_vector = torch.column_stack((torch.tensor([1, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, 1, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, 1, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, 1, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, 1, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 1,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,1,-np.inf,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,1,-np.inf,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,1,-np.inf,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,1,-np.inf]),
                                 torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf,-np.inf,-np.inf,-np.inf,-np.inf,1])
                                   ))
    
    out = loss_model._reduce_logits(input_vector)
    assert torch.equal(out, torch.tensor([[1, -np.inf, -np.inf, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, 1, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, -np.inf, 1],
                                          [-np.inf, -np.inf, 1, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, 1, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, -np.inf, 1],
                                          [-np.inf, 1, -np.inf, -np.inf, -np.inf],
                                          [-np.inf, 1, -np.inf, -np.inf, -np.inf],
                                          [-np.inf, 1, -np.inf, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, 1, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, 1, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, -np.inf, 1],
                                          [-np.inf, -np.inf, 1, -np.inf, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, 1, -np.inf],
                                          [-np.inf, -np.inf, -np.inf, -np.inf, 1],
                                          ])), "Reduced logit vector did not match expected"

def test_predict_codon():
    loss_model = loss.TiberiusLoss(use_f1=True, f1_factor=2.0)
    #first, test the one-hots
    assert loss_model._predict_codon(torch.tensor([1, -np.inf, -np.inf, -np.inf, -np.inf])) == False
    assert loss_model._predict_codon(torch.tensor([-np.inf, 1, -np.inf, -np.inf, -np.inf])) == False
    assert loss_model._predict_codon(torch.tensor([-np.inf, -np.inf, 1, -np.inf, -np.inf])) == True
    assert loss_model._predict_codon(torch.tensor([-np.inf, -np.inf, -np.inf, 1, -np.inf])) == True
    assert loss_model._predict_codon(torch.tensor([-np.inf, -np.inf, -np.inf, -np.inf, 1])) == True

    #now, some combos
    assert loss_model._predict_codon(torch.tensor([1, 1, -np.inf, 1, -np.inf])) == False
    assert loss_model._predict_codon(torch.tensor([2, -np.inf, -np.inf, -np.inf, 1])) == False
    assert loss_model._predict_codon(torch.tensor([1, -np.inf, 1, 1, -np.inf])) == True
