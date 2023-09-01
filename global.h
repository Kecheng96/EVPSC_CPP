#ifndef GLOBALS_H
#define GLOBALS_H
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <ctime>
#include <Eigen/Dense>

extern double temp_K;
extern std::ofstream outfile;
class Logger;
extern Logger logger; // Global logger object, include global.h to use it

void update_progress(double progress_f);

// Logger class for logging messages to the console and a log file
class Logger {
public:
    Logger() : outfile("EVPSC_log.txt") {}

    // Define the message hierarchy levels inside the Logger class
    enum class MessageLevel {
        DEBUG,
        NOTICE,
        WARN,
        INFO,
        ERROR,
    };

    // Log a message at a specific message level
    void log(MessageLevel level, const std::string& message) {
        // Get the current time
        std::time_t t = std::time(nullptr);
        std::string timestamp = std::asctime(std::localtime(&t));
        timestamp.pop_back();  // Remove the newline character from the end

        // Output the message to the console
        if (console_level_ <= level) {
            std::cout << "[" << timestamp << "] ";
            printPrefix(level);
            std::cout << message << std::endl;
        }

        // Output the message to the log file
        if (file_level_ <= level) {
            outfile << "[" << timestamp << "] ";
            printPrefix(level, outfile);
            outfile << message << std::endl;
        }
    }

    void log(MessageLevel level, const Eigen::MatrixXd& message) {
        // Get the current time
        std::time_t t = std::time(nullptr);
        std::string timestamp = std::asctime(std::localtime(&t));
        timestamp.pop_back();  // Remove the newline character from the end
        Eigen::IOFormat LogFmt(4, 0, ", ", "\n", "["+timestamp+"] "+getPrefix(level)+ "[", "]");
        // Output the message to the console
        if (console_level_ <= level) {
            std::cout << message.format(LogFmt) << std::endl;
        }
        // Output the message to the log file
        if (file_level_ <= level) {
            outfile << message.format(LogFmt) << std::endl;
        }
    }

    // Overloaded functions for different message levels
    void debug(const std::string& message) { log(MessageLevel::DEBUG, message); }
    void info(const std::string& message) { log(MessageLevel::INFO, message); }
    void notice(const std::string& message) { log(MessageLevel::NOTICE, message); }
    void warn(const std::string& message) { log(MessageLevel::WARN, message); }
    void error(const std::string& message) { log(MessageLevel::ERROR, message); }
    void debug(const Eigen::MatrixXd& message) { log(MessageLevel::DEBUG, message); }
    void info(const Eigen::MatrixXd& message) { log(MessageLevel::INFO, message); }
    void notice(const Eigen::MatrixXd& message) { log(MessageLevel::NOTICE, message); }
    void warn(const Eigen::MatrixXd& message) { log(MessageLevel::WARN, message); }
    void error(const Eigen::MatrixXd& message) { log(MessageLevel::ERROR, message); }
    // Set the message level filter for console output
    void setConsoleLevel(MessageLevel level) { console_level_ = level; }

    // Set the message level filter for file output
    void setFileLevel(MessageLevel level) { file_level_ = level; }

private:
    std::ofstream outfile;
    MessageLevel console_level_ = MessageLevel::INFO;
    MessageLevel file_level_ = MessageLevel::DEBUG;

    // Helper function to print the message prefix based on the message level
    void printPrefix(MessageLevel level, std::ostream& stream = std::cout) {
        switch (level) {
            case MessageLevel::DEBUG:
                stream << "[ DEBUG] ";
                break;
            case MessageLevel::INFO:
                stream << "[ INFO ] ";
                break;
            case MessageLevel::NOTICE:
                stream << "[NOTICE] ";
                break;
            case MessageLevel::WARN:
                stream << "[ WARN ] ";
                break;
            case MessageLevel::ERROR:
                stream << "[ ERROR] ";
                break;
            default:
                break;
        }
    }
    std::string getPrefix(MessageLevel level) {
        switch (level) {
            case MessageLevel::DEBUG:
                return "[ DEBUG] ";
            case MessageLevel::INFO:
                return "[ INFO ] ";
            case MessageLevel::NOTICE:
                return "[NOTICE] ";
            case MessageLevel::WARN:
                return "[ WARN ] ";
            case MessageLevel::ERROR:
                return "[ ERROR] ";
            default:
                return "";
        }
    }
};


#endif

