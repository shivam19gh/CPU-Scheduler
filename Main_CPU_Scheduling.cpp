#include <iostream>
#include <queue>
#include <utility>
#include <string>
#include <vector>
#include <algorithm> // for all_of, count_if
#include <numeric>  // for accumulate
#include <cmath>
#include <limits>
#include <map>
#include <fstream>
using namespace std;

struct Process{
    string pid; // prcess id
    int bt; // burst time
    int at; // arrival time
    int priority; // priority of the process
    int wt; // wait time
    int tat; // turn around time
    Process(){
      pid = ""; bt = 0; at =0; priority =0; wt =0; tat =0;
    }
};

bool priorityVary(vector<Process> &processes){
   int pi = processes[0].priority;
   for(auto &it: processes){
    if(it.priority!=pi) return true;
   }
return false;
}

string determineSchedulingAlgorithm(vector<Process> &tasks, int high_threshold, int low_threshold, double &burst_time_variance, double &arrival_time_variance, double &average_burst_time) {
    // Compute total burst time for all processes
    double total_burst_time = accumulate(tasks.begin(), tasks.end(), 0.0, [](double sum, Process &p) { return sum + p.bt; });
    average_burst_time = total_burst_time / tasks.size();

    // Check if priority varies among processes
    bool has_priority_variation = priorityVary(tasks);

    // Calculate variance of burst times
    burst_time_variance = 0;
    for (const auto &task : tasks) {
        burst_time_variance += pow(task.bt - average_burst_time, 2);
    }
    burst_time_variance /= tasks.size();

    // Compute total and average arrival times
    double total_arrival_time = accumulate(tasks.begin(), tasks.end(), 0.0, [](double sum, const Process &p) { return sum + p.at; });
    double average_arrival_time = total_arrival_time / tasks.size();

    // Calculate variance of arrival times
    arrival_time_variance = 0;
    for (const auto &task : tasks) {
        arrival_time_variance += pow(task.at - average_arrival_time, 2);
    }
    arrival_time_variance /= tasks.size();

    // Determine scheduling algorithm based on variance and priority
    if (has_priority_variation) {
        return "priorityBasedScheduling";
    } else if (arrival_time_variance > high_threshold) {
        if (burst_time_variance > high_threshold) {
            return "shortestRemainingTimeFirst";
        } else {
            return "firstComeFirstServe";
        }
    } else if (burst_time_variance > high_threshold) {
        return "shortestJobFirst";
    } else if (burst_time_variance < low_threshold && average_burst_time > high_threshold) {
        return "longestJobFirst";
    } else {
        return "roundRobinScheduling";
    }
}

// Scheduling algorithms
void firstComeFirstServe(vector<Process> &tasks, double burstTimeVariance, double arrivalTimeVariance, double avgBurstTime) {
    // First Come First Serve (FCFS) is a non-preemptive scheduling algorithm, processes are executed in the order of their arrival times.
    int currentTime = 0;
    double totalTurnaroundTime = 0, totalWaitTime = 0, totalResponseTime = 0;

    for (int i = 0; i < tasks.size(); i++) {
        if (tasks[i].at > currentTime) {
            currentTime = tasks[i].at;
        }

        totalTurnaroundTime += static_cast<double>(currentTime + tasks[i].bt - tasks[i].at);
        totalWaitTime += static_cast<double>(currentTime - tasks[i].at);
        totalResponseTime += static_cast<double>(currentTime - tasks[i].at);

        currentTime += tasks[i].bt;
    }

    int numProcesses = tasks.size();
    double avgTurnaroundTime = totalTurnaroundTime / static_cast<double>(numProcesses);
    double avgWaitTime = totalWaitTime / static_cast<double>(numProcesses);
    double avgResponseTime = totalResponseTime / static_cast<double>(numProcesses);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
      outputFile << "The FCFS algorithm is ideal for scenarios where there is significant variability in arrival times." << endl;
      outputFile << "Arrival time variance: " << to_string(arrivalTimeVariance) << endl;
      outputFile << "Burst time variance: " << to_string(burstTimeVariance) << endl;
      outputFile << "The average turnaround time for the given processes is: " << to_string(avgTurnaroundTime) << endl;
      outputFile << "The average waiting time for the given processes is: " << to_string(avgWaitTime) << endl;
      outputFile << "The average response time for the given processes is: " << to_string(avgResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process Execution Order: " << endl;
        for (int i = 0; i < tasks.size(); i++) {
            outputFile << tasks[i].pid << " ";
        }
        outputFile.close();
    }
}


void shortestJobFirst(vector<Process> &tasks, double burstTimeVariance, double arrivalTimeVariance, double avgBurstTime) {
    // SJF is a generally non-preemptive scheduling algorithm
    vector<string> executionOrder;
    int currentTime = tasks[0].at, taskIndex = 0, totalTasks = tasks.size();
    double totalTurnaroundTime = 0, totalWaitTime = 0, totalResponseTime = 0;
    priority_queue<pair<int, Process*>, vector<pair<int, Process*>>, greater<pair<int, Process*>>> readyQueue;

    Process* currentTask = &tasks[taskIndex];
    readyQueue.push({tasks[taskIndex].bt, currentTask});
    taskIndex++;

    while (!readyQueue.empty()) {

        while (taskIndex < totalTasks && tasks[taskIndex].at <= currentTime) {
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }

        pair<int, Process*> nextTask = readyQueue.top();  
        readyQueue.pop();
        executionOrder.push_back(nextTask.second->pid);
        currentTime += nextTask.first;
        nextTask.second->tat = currentTime - nextTask.second->at;
        nextTask.second->wt = nextTask.second->tat - nextTask.second->bt;
        totalTurnaroundTime += nextTask.second->tat;
        totalWaitTime += nextTask.second->wt;
        totalResponseTime += currentTime - nextTask.first - nextTask.second->at;

        if (readyQueue.empty() && taskIndex < totalTasks) {
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }
    }

    double avgTurnaroundTime = totalTurnaroundTime / static_cast<double>(totalTasks);
    double avgWaitTime = totalWaitTime / static_cast<double>(totalTasks);
    double avgResponseTime = totalResponseTime / static_cast<double>(totalTasks);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
        outputFile << "SJF excels in reducing the average wait time when there is significant variability in burst durations." << endl;
        outputFile << "Variance in arrival times: " << to_string(arrivalTimeVariance) << endl;
        outputFile << "Variance in burst times: " << to_string(burstTimeVariance) << endl;
        outputFile << "Average turnaround time for the provided dataset: " << to_string(avgTurnaroundTime) << endl;
        outputFile << "Average waiting time for the provided dataset: " << to_string(avgWaitTime) << endl;
        outputFile << "Average response time for the provided dataset: " << to_string(avgResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process Execution Order: " << endl;
        for (const auto& pid : executionOrder) {
            outputFile << pid << " ";
        }
        outputFile.close();
    }
}

void longestJobFirst(vector<Process> &tasks, double burstTimeVariance, double arrivalTimeVariance, double avgBurstTime) {
    vector<string> executionOrder;
    int currentTime = tasks[0].at, taskIndex = 0, totalTasks = tasks.size();
    double totalTurnaroundTime = 0, totalWaitTime = 0, totalResponseTime = 0;
    priority_queue<pair<int, Process*>> readyQueue;

    Process* currentTask = &tasks[taskIndex];
    readyQueue.push({tasks[taskIndex].bt, currentTask});
    taskIndex++;

    while (!readyQueue.empty()) {
        // Add tasks to the ready queue whose arrival times are less than or equal to the current CPU time
        while (taskIndex < totalTasks && tasks[taskIndex].at <= currentTime) {
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }

        pair<int, Process*> nextTask = readyQueue.top();  // Next task to process
        readyQueue.pop();
        executionOrder.push_back(nextTask.second->pid);
        currentTime += nextTask.first;
        nextTask.second->tat = currentTime - nextTask.second->at;
        nextTask.second->wt = nextTask.second->tat - nextTask.second->bt;
        totalTurnaroundTime += nextTask.second->tat;
        totalWaitTime += nextTask.second->wt;
        totalResponseTime += currentTime - nextTask.first - nextTask.second->at;

        if (readyQueue.empty() && taskIndex < totalTasks) {
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }
    }

    double avgTurnaroundTime = totalTurnaroundTime / static_cast<double>(totalTasks);
    double avgWaitTime = totalWaitTime / static_cast<double>(totalTasks);
    double avgResponseTime = totalResponseTime / static_cast<double>(totalTasks);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
        outputFile << "LJFS gives preference to longer tasks when the variability in burst times is minimal." << endl;
        outputFile << "Arrival time variance: " << to_string(arrivalTimeVariance) << endl;
        outputFile << "Burst time variance: " << to_string(burstTimeVariance) << endl;
        outputFile << "Average turnaround time for the provided dataset: " << to_string(avgTurnaroundTime) << endl;
        outputFile << "Average waiting time for the provided dataset: " << to_string(avgWaitTime) << endl;
        outputFile << "Average response time for the provided dataset: " << to_string(avgResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process execution Order : " << endl;
        for (const auto& pid : executionOrder) {
            outputFile << pid << " ";
        }
        outputFile.close();
    }
}

void shortestRemainingTimeFirst(vector<Process> &tasks, double burstTimeVariance, double arrivalTimeVariance, double avgBurstTime) {
    vector<string> executionOrder;
    int currentTime = tasks[0].at, taskIndex = 0, totalTasks = tasks.size();
    double totalTurnaroundTime = 0, totalWaitTime = 0, totalResponseTime = 0;
    priority_queue<pair<int, Process*>, vector<pair<int, Process*>>, greater<pair<int, Process*>>> readyQueue;

    Process* currentTask = &tasks[taskIndex];
    readyQueue.push({tasks[taskIndex].bt, currentTask});
    taskIndex++;

    while (!readyQueue.empty()) {
        // Add tasks to the ready queue whose arrival times are less than or equal to the current CPU time
        while (taskIndex < totalTasks && tasks[taskIndex].at <= currentTime) {
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }

        pair<int, Process*> nextTask = readyQueue.top();  // Next task to process
        readyQueue.pop();
        if (executionOrder.empty() || executionOrder.back() != nextTask.second->pid) {
            executionOrder.push_back(nextTask.second->pid);
        }

        if (nextTask.first == nextTask.second->bt) {
            totalResponseTime += (currentTime - nextTask.second->at);
        }
        currentTime += 1;
        nextTask.first -= 1;
        if (nextTask.first == 0) {
            nextTask.second->tat = currentTime - nextTask.second->at;
            nextTask.second->wt = nextTask.second->tat - nextTask.second->bt;
            totalTurnaroundTime += nextTask.second->tat;
            totalWaitTime += nextTask.second->wt;
        } else {
            readyQueue.push(nextTask);
        }

        if (readyQueue.empty() && taskIndex < totalTasks) {
            currentTime = tasks[taskIndex].at;
            currentTask = &tasks[taskIndex];
            readyQueue.push({tasks[taskIndex].bt, currentTask});
            taskIndex++;
        }
    }

    double avgTurnaroundTime = totalTurnaroundTime / static_cast<double>(totalTasks);
    double avgWaitTime = totalWaitTime / static_cast<double>(totalTasks);
    double avgResponseTime = totalResponseTime / static_cast<double>(totalTasks);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
        outputFile << "Due to high variability in both arrival and burst times, SRTF is selected to efficiently manage dynamic changes." << endl;
        outputFile << "Variance in arrival times: " << to_string(arrivalTimeVariance) << endl;
        outputFile << "Variance in burst times: " << to_string(burstTimeVariance) << endl;
        outputFile << "Average turnaround time for the provided dataset: " << to_string(avgTurnaroundTime) << endl;
        outputFile << "Average waiting time for the provided dataset: " << to_string(avgWaitTime) << endl;
        outputFile << "Average response time for the provided dataset: " << to_string(avgResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process execution Order : " << endl;
        for (const auto& pid : executionOrder) {
            outputFile << pid << " ";
        }
        outputFile.close();
    }
}

void roundRobinScheduling(vector<Process> &processes, double burstTimeVariance, double arrivalTimeVariance, double averageBurstTime) {
    vector<string> executionOrder;
    int currentTime = processes[0].at, i = 0, processCount = processes.size();
    double totalTurnaroundTime = 0, totalWaitTime = 0, totalResponseTime = 0;
    queue<pair<int, Process*>> readyQueue;
    Process* currentProcess = &processes[0];
    readyQueue.push({currentProcess->bt, currentProcess});
    i++;

    while (!readyQueue.empty()) {
        pair<int, Process*> currentTask = readyQueue.front();
        readyQueue.pop();

        if (currentTask.first == currentTask.second->bt) {
            totalResponseTime += (currentTime - currentTask.second->at);
        }

        if (executionOrder.empty() || executionOrder.back() != currentTask.second->pid) {
            executionOrder.push_back(currentTask.second->pid);
        }

        currentTime += 1;
        currentTask.first -= 1;

        while (i < processCount && processes[i].at <= currentTime) {
            currentProcess = &processes[i];
            readyQueue.push({currentProcess->bt, currentProcess});
            i++;
        }

        if (currentTask.first == 0) {
            totalTurnaroundTime += currentTime - currentTask.second->at;
            totalWaitTime += (currentTime - currentTask.second->at - currentTask.second->bt);
        } else {
            readyQueue.push(currentTask);
        }

        if (readyQueue.empty() && i < processCount) {
            currentProcess = &processes[i];
            readyQueue.push({currentProcess->bt, currentProcess});
            currentTime = currentProcess->at;
            i++;
        }
    }

    totalTurnaroundTime /= double(processCount);
    totalWaitTime /= double(processCount);
    totalResponseTime /= double(processCount);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
        outputFile << "Round Robin is chosen for equitable processing and responsiveness in time-sharing environments with minimal variation in arrival and burst times." << endl;
        outputFile << "Variance of arrival times: " << to_string(arrivalTimeVariance) << endl;
        outputFile << "Variance of burst times: " << to_string(burstTimeVariance) << endl;
        outputFile << "Average turnaround time: " << to_string(totalTurnaroundTime) << endl;
        outputFile << "Average wait time: " << to_string(totalWaitTime) << endl;
        outputFile << "Average response time: " << to_string(totalResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process execution Order : " << endl;
        for (const auto &pid : executionOrder) outputFile << pid << " ";
        outputFile.close();
    }
}

void priorityBasedScheduling(vector<Process> &processes, double burstTimeVariance, double arrivalTimeVariance, double averageBurstTime) {
    vector<string> executionOrder;
    int currentTime = 0, i = 0, processCount = processes.size();
    double totalTurnaroundTime = 0, totalResponseTime = 0, totalWaitTime = 0;
    priority_queue<pair<int, pair<int, Process*>>> priorityQueue;
    currentTime = processes[i].at;
    Process* currentProcess;

    while (i < processCount && processes[i].at <= currentTime) {
        currentProcess = &processes[i];
        priorityQueue.push({processes[i].priority, {currentProcess->bt, currentProcess}});
        i++;
    }

    while (!priorityQueue.empty()) {
        pair<int, pair<int, Process*>> currentTask = priorityQueue.top();
        priorityQueue.pop();

        if (currentTask.second.first == currentTask.second.second->bt) {
            totalResponseTime += (currentTime - currentTask.second.second->at);
        }

        if (executionOrder.empty() || executionOrder.back() != currentTask.second.second->pid) {
            executionOrder.push_back(currentTask.second.second->pid);
        }

        currentTime += 1;
        currentTask.second.first -= 1;

        if (currentTask.second.first == 0) {
            totalTurnaroundTime += currentTime - currentTask.second.second->at;
            totalWaitTime += (currentTime - currentTask.second.second->at - currentTask.second.second->bt);
        } else {
            priorityQueue.push(currentTask);
        }

        while (i < processCount && processes[i].at <= currentTime) {
            currentProcess = &processes[i];
            priorityQueue.push({currentProcess->priority, {currentProcess->bt, currentProcess}});
            i++;
        }

        if (i < processCount && priorityQueue.empty()) {
            currentTime = processes[i].at;
            while (i < processCount && processes[i].at <= currentTime) {
                currentProcess = &processes[i];
                priorityQueue.push({currentProcess->priority, {currentProcess->bt, currentProcess}});
                i++;
            }
        }
    }

    totalTurnaroundTime /= double(processCount);
    totalWaitTime /= double(processCount);
    totalResponseTime /= double(processCount);

    fstream outputFile;
    outputFile.open("output.txt", ios::out);
    if (outputFile.is_open()) {
        outputFile << "Priority Scheduling has been selected as the scheduling algorithm." << endl;
        outputFile << "Average turnaround time: " << to_string(totalTurnaroundTime) << endl;
        outputFile << "Average wait time: " << to_string(totalWaitTime) << endl;
        outputFile << "Average response time: " << to_string(totalResponseTime) << endl;
        outputFile << endl;
        outputFile << "Process execution Order : " << endl;
        for (const auto &pid : executionOrder) outputFile << pid << " ";
        outputFile.close();
    }
}


bool isValidProcess(Process &p) {
    if (p.bt == 0 || p.pid.empty())
        return false;
    return true;
}

void getData(vector<Process> &processes){
     fstream inputFile;
     inputFile.open("input.txt",ios::in);
     while(inputFile.is_open()){
        string pdata;
		while(getline(inputFile,pdata)){
            Process p;
			string data=""; int nd=0,i=0,n = pdata.length();
        while(nd<4&&i<n){
          if(pdata[i]==' '){
				   i++;
				   continue;
			    }
          else{
				  nd++;
				  while(pdata[i]!=' '&&i<n){
					data +=pdata[i];
					i++;
				   } 
                // if nd =1 -> pid, nd=2 -> arrival time , nd =3 -> burst time, nd = 4 -> priority
                if(nd==1) p.pid = data;
                else if(nd ==2) p.at = stoi(data);
                else if(nd == 3) p.bt = stoi(data);
                else p.priority = stoi(data);
               data = "";
			   }
			}
      if(isValidProcess(p))
      processes.push_back(p);
		}
      inputFile.close();
     }
   sort(processes.begin(),processes.end(),[](Process &a, Process &b){return a.at<=b.at;});
}

int main() {
    int variance_threshold_high = 10, variance_threshold_low = 2;
    double vbt = 0, vat = 0, abt = 0;
    vector<Process> processes;
    
    // Function to retrieve process data from input
    getData(processes);
    
    // Determining the scheduling algorithm to apply
    string algorithmToApply = determineSchedulingAlgorithm(processes, variance_threshold_high, variance_threshold_low, vbt, vat, abt);
    
    // Applying the selected scheduling algorithm
    if(algorithmToApply.compare("firstComeFirstServe")==0) firstComeFirstServe(processes, vbt, vat, abt);
    else if(algorithmToApply.compare("shortestJobFirst")==0) shortestJobFirst(processes, vbt, vat, abt);
    else if(algorithmToApply.compare("longestJobFirst")==0) longestJobFirst(processes, vbt, vat, abt);
    else if(algorithmToApply.compare("shortestRemainingTimeFirst")==0) shortestRemainingTimeFirst(processes, vbt, vat, abt);
    else if(algorithmToApply.compare("roundRobinScheduling")==0) roundRobinScheduling(processes, vbt, vat, abt);
    else priorityBasedScheduling(processes, vbt, vat, abt);
    return 0;
}
