#define CROW_MAIN
#include "crow_all.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

// Function to run the scheduling algorithm
void runSchedulingAlgorithm() {
    system("Main_CPU_Scheduling.exe");
}

int main() {
    crow::SimpleApp app;

    // Route for the main page to upload the file
    CROW_ROUTE(app, "/").methods(crow::HTTPMethod::GET)([]() {
        std::ostringstream html;
        html << "<html><body>"
             << "<h1>Upload CPU Scheduling Input File</h1>"
             << "<form action='/upload' method='post' enctype='multipart/form-data'>"
             << "<input type='file' name='file'>"
             << "<input type='submit' value='Upload'>"
             << "</form></body></html>";
        return crow::response(html.str());
    });

    // Route to handle the file upload
    CROW_ROUTE(app, "/upload").methods(crow::HTTPMethod::POST)([](const crow::request& req) {
        auto file = req.get_part("file");

        if (!file) {
            return crow::response(400, "No file uploaded");
        }

        // Save the uploaded file to 'input.txt'
        std::ofstream ofs("input.txt", std::ios::binary);
        if (ofs.is_open()) {
            ofs.write(file->body.data(), file->body.size());
            ofs.close();
        } else {
            return crow::response(500, "Failed to save the uploaded file");
        }

        // Run the scheduling algorithm
        runSchedulingAlgorithm();

        // Read the output from 'output.txt'
        std::ifstream ifs("output.txt");
        if (!ifs.is_open()) {
            return crow::response(500, "Failed to read the output file");
        }

        std::stringstream buffer;
        buffer << ifs.rdbuf();
        ifs.close();

        return crow::response(buffer.str());
    });

    // Start the Crow server on port 8080
    app.port(8080).multithreaded().run();
}