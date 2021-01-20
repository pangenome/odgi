#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <thread>
#include "atomic_queue.h"

namespace odgi {
namespace algorithms {

class bed_records_class {

private:
    // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
    struct bed_record_t {
        std::string chrom;
        uint64_t chromStart;
        uint64_t chromEnd;
        double path_layout_dist;
        uint64_t path_nuc_dist;
    };

    std::thread writer_thread;
    atomic_queue::AtomicQueue2<bed_record_t*, 2 << 16> bed_record_queue;
    std::atomic<bool> work_todo;

public:

    bed_records_class(void) {
        work_todo.store(false);
    }

    ~bed_records_class(void) {
        close_writer();
    }

    void writer_func(void) {
        /*
        std::ofstream writer(filename.c_str(), std::ios::binary | std::ios::trunc);
        if (writer.fail()) {
            throw std::ios_base::failure(std::strerror(errno));
        }
        */
        auto* bed_record = new bed_record_t();
        while (work_todo.load() || !bed_record_queue.was_empty()) {
            if (bed_record_queue.try_pop(bed_record)) {
                do {
                    // writer.write((char*)&ival, sizeof(Interval));
                    // BED files are 0-based http://genome.ucsc.edu/FAQ/FAQformat#format1
                    std::cout << bed_record->chrom << "\t" // chrom
                              << bed_record->chromStart << "\t" // chromStart
                              << bed_record->chromEnd << "\t" // chromEnd
                              << bed_record->path_layout_dist << "\t"
                              << bed_record->path_nuc_dist
                              << std::endl;
                } while (bed_record_queue.try_pop(bed_record));
            } else {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
        }
        // writer.close();
    }

    // open backing file and start writer_thread
    void open_writer(void) {
        std::cerr << "PENEN" << std::endl;
        if (!work_todo.load()) {
            work_todo.store(true);
            writer_thread = std::thread(&bed_records_class::writer_func, this);
        }
    }

    void close_writer(void) {
        if (work_todo.load()) {
            work_todo.store(false);
            while (!bed_record_queue.was_empty()) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            if (writer_thread.joinable()) {
                writer_thread.join();
            }
        }
    }

    /// write into our write buffer
    /// open_writer() must be called first to set up our buffer and writer
	void append(const std::string &chrom, const uint64_t &chromStart, const uint64_t &chromEnd,
          const double &path_layout_dist, const uint64_t &path_nuc_dist) {
        std::cerr << "sdfa" << std::endl;
        bed_record_t bed_record = {chrom, chromStart, chromEnd, path_layout_dist, path_nuc_dist};
        bed_record_queue.push(&bed_record);
    }
};

}

}