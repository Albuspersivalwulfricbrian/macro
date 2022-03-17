//
// Created by Sultan Musin on 26.10.2021.
//

#ifndef ENT_ROOT_TREE_STRUCTURES_H
#define ENT_ROOT_TREE_STRUCTURES_H


struct EnergyCollection {
    Float_t scat0;
    Float_t scat1;
    Float_t det0;
    Float_t det1;
    Float_t intermediate;

    Bool_t operator==(const EnergyCollection &sec) const {
        return scat0 == sec.scat0 && scat1 == sec.scat1 && det0 == sec.det0 &&
               det1 == sec.det1 && intermediate == sec.intermediate;
    }

    Bool_t operator!=(const EnergyCollection &sec) const {
        return !(operator==(sec));
    }

    Bool_t operator<=(const EnergyCollection &sec) const {
        return scat0 <= sec.scat0 && scat1 <= sec.scat1 && det0 <= sec.det0
               && det1 <= sec.det1 && intermediate <= sec.intermediate;
    };

    Bool_t operator<(const EnergyCollection &sec) const {
        return scat0 < sec.scat0 && scat1 < sec.scat1 && det0 < sec.det0 &&
               det1 < sec.det1 && intermediate < sec.intermediate;;
    }

    Bool_t operator>=(const EnergyCollection &sec) const {
        return !(operator<(sec));
    }

    Bool_t operator>(const EnergyCollection &sec) const {
        return !(operator<=(sec));
    }
};

class TreeStructures {
public:
    explicit TreeStructures() {};

    ~TreeStructures() {};
};

class InitialTreeStructure : public TreeStructures {
public:
    InitialTreeStructure() : TreeStructures() {};

    Short_t Get_wf_size() { return wf_size; }

    Short_t *Get_wf() { return wf; };

private:
    Short_t wf_size;
    Short_t wf[2048];
};

class CalibratedTreeStructure : public TreeStructures {
public:
    CalibratedTreeStructure() : TreeStructures(),
                                integral_in_gate(0), peak_pos(0), amp(0) {};

    CalibratedTreeStructure(Int_t integral, Int_t peak, Int_t am) :
            integral_in_gate(integral), peak_pos(peak), amp(am) {};

    Int_t Get_integral_in_gate() { return integral_in_gate; };

    Int_t Get_peak_pos() { return peak_pos; };

    Int_t Get_amp() { return amp; };

    void CopyMembers(CalibratedTreeStructure const &ne) {
        this->integral_in_gate = ne.integral_in_gate;
        this->peak_pos = ne.peak_pos;
        this->amp = ne.amp;
    };

private:
    Int_t integral_in_gate;
    Int_t peak_pos;
    Int_t amp;

};

class MiniTree : public TreeStructures {
public:
    MiniTree() : TreeStructures() {};

    MiniTree(Float_t Inter, Float_t Scat0, Float_t Scat1, Float_t Det0, Float_t Det1, Short_t DetN0, Short_t DetN1) :
            EdepIntermediate(Inter), EdepScat0(Scat0), EdepScat1(Scat1),
            EdepDet0(Det0), EdepDet1(Det1), DetNum0(DetN0), DetNum1(DetN1) {};

    Float_t GetEdepIntermediate() { return EdepIntermediate; };

    Float_t GetEdepScat0() { return EdepScat0; };

    Float_t GetEdepScat1() { return EdepScat1; };

    Float_t GetEdepDet0() { return EdepDet0; };

    Float_t GetEdepDet1() { return EdepDet1; };

    Short_t GetDetNum0() { return DetNum0; };

    Short_t GetDetNum1() { return DetNum1; };

    void CopyMembers(MiniTree const &ne) {
        this->EdepIntermediate = ne.EdepIntermediate;
        this->EdepScat0 = ne.EdepScat0;
        this->EdepScat1 = ne.EdepScat1;
        this->EdepDet0 = ne.EdepDet0;
        this->EdepDet1 = ne.EdepDet1;
        this->DetNum0 = ne.DetNum0;
        this->DetNum1 = ne.DetNum1;
    };

    Short_t DistanceBetweenNaI() {
        Short_t angle;
        angle = DetNum1 - 16 - DetNum0;
        if (angle < 0) angle += 16;
        return angle;
    }

    std::array<Float_t, 5> GetEdepArray() {
        return {EdepScat0, EdepScat1, EdepDet0, EdepDet1, EdepIntermediate};
    };

    EnergyCollection GetEnergyCollection() {
        return {EdepScat0, EdepScat1, EdepDet0, EdepDet1, EdepIntermediate};
    }

public:
    Float_t EdepIntermediate;
    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;
    Short_t DetNum0;
    Short_t DetNum1;
};


#endif //ENT_ROOT_TREE_STRUCTURES_H
