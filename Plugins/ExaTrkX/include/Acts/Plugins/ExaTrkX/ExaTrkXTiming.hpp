#pragma once

#include <vector>
#include <numeric>
#include <chrono>
#include <ctime>
#include <stdio.h>
#include <string>

struct ExaTrkXTime {
    float embedding = 0.0;
    float building = 0.0;
    float filtering = 0.0;
    float gnn = 0.0;
    float labeling = 0.0;
    float total = 0.0;
    void summary() {
        printf("1) embedding:  %.4f\n", embedding);
        printf("2) building:   %.4f\n", building);
        printf("3) filtering:  %.4f\n", filtering);
        printf("4) gnn:        %.4f\n", gnn);
        printf("5) labeling:   %.4f\n", labeling);
        printf("6) total:      %.4f\n", total);
    }
};

struct ExaTrkXTimeList {
    std::vector<float> embedding;
    std::vector<float> building;
    std::vector<float> filtering;
    std::vector<float> gnn;
    std::vector<float> labeling;
    std::vector<float> total;

    void add(const ExaTrkXTime& time) {
        embedding.push_back(time.embedding);
        building.push_back(time.building);
        filtering.push_back(time.filtering);
        gnn.push_back(time.gnn);
        labeling.push_back(time.labeling);
        total.push_back(time.total);
    }

    ExaTrkXTime get(int evtid) {
        if (evtid >= embedding.size()) {
            printf("Error: event id %d is out of range.\n", evtid);
            return ExaTrkXTime();
        }
        ExaTrkXTime timing {
            embedding[evtid],
            building[evtid],
            filtering[evtid],
            gnn[evtid],
            labeling[evtid],
            total[evtid]
        };
        return timing;
    }

    void summary(int start=0) {
        size_t num = embedding.size();
        if (num <= start) {
            printf("Not enough data. %ld total and %d skipped\n", num, start);
            return;
        }
        num -= start;
        float tot_embedding = 0;
        float tot_building = 0;
        float tot_filtering = 0;
        float tot_gnn = 0;
        float tot_labeling = 0;
        float tot_total = 0;
        if (num > 0){
            tot_embedding = std::accumulate(embedding.begin()+start, embedding.end(), 0.0f);
            tot_building = std::accumulate(building.begin()+start, building.end(), 0.0f);
            tot_filtering = std::accumulate(filtering.begin()+start, filtering.end(), 0.0f);
            tot_gnn = std::accumulate(gnn.begin()+start, gnn.end(), 0.0f);
            tot_labeling = std::accumulate(labeling.begin()+start, labeling.end(), 0.0f);
            tot_total = std::accumulate(total.begin()+start, total.end(), 0.0f);
        }

        printf("1) embedding: %.4f\n", tot_embedding / num);
        printf("2) building:  %.4f\n", tot_building / num);
        printf("3) filtering: %.4f\n", tot_filtering / num);
        printf("4) gnn:       %.4f\n", tot_gnn / num);
        printf("5) labeling:  %.4f\n", tot_labeling / num);
        printf("6) total:     %.4f\n", tot_total / num);
    }

    void summaryOneEvent(int evtid){
        get(evtid).summary();
    }
    
    int numEvts(){
        return (int) embedding.size();
    }

    void save(const std::string& filename) {
        FILE* fp = fopen(filename.c_str(), "w");
        if (fp == NULL) {
            printf("Error: cannot open file %s\n", filename.c_str());
            return;
        }
        fprintf(fp, "embedding,building,filtering,gnn,labeling,total\n");
        for (int i = 0; i < embedding.size(); i++) {
            fprintf(fp, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                embedding[i], building[i], filtering[i], gnn[i], labeling[i], total[i]);
        }
        fclose(fp);
    }


};

class ExaTrkXTimer
{
public:
    void start() { 
        m_start = std::chrono::high_resolution_clock::now(); m_running = true;
    }
    void stop() {
        m_end = std::chrono::high_resolution_clock::now(); m_running = false;
    }
    double stopAndGetElapsedTime() {
        stop();
        return elapsedSeconds();
    }
    double elapsed() {
        std::chrono::time_point<std::chrono::high_resolution_clock> end;
        if (m_running) {
            end = std::chrono::high_resolution_clock::now();
        } else { end = m_end; }
        return std::chrono::duration<double, std::milli>(end - m_start).count();
    }
    double elapsedSeconds() {
        return elapsed() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
    bool m_running = false;
};