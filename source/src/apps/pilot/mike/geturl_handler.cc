// Copyright (c) 2012 The Chromium Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include <stdio.h>
#include <stdlib.h>
#include "ppapi/c/pp_errors.h"
#include "ppapi/c/ppb_instance.h"
#include "ppapi/cpp/module.h"
#include "ppapi/cpp/var.h"
#include "ppapi/cpp/url_response_info.h"

#include "geturl_handler.h"
#include <iostream>
#include <stringstream>

#ifdef WIN32
#undef min
#undef max
#undef PostMessage
#endif

GetURLHandler* GetURLHandler::Create(pp::Instance* instance,
                                     const std::string& url,
                                     utility::Inline_File_Provider_Hook* file_provider_hook ) {
  return new GetURLHandler(instance, url, file_provider_hook);
}

GetURLHandler::GetURLHandler(pp::Instance* instance,
                             const std::string& url,
                             utility::Inline_File_Provider_Hook* file_provider_hook  
                            )
    : instance_(instance),
      url_(url),
      url_request_(instance),
      url_loader_(instance),
      buffer_(new char[READ_BUFFER_SIZE]),
      cc_factory_(this),
      file_provider_hook_(file_provider_hook) {
  
  url_request_.SetURL(url);
  url_request_.SetMethod("GET");
  url_request_.SetRecordDownloadProgress(true);
}

GetURLHandler::~GetURLHandler() {
  delete [] buffer_;
  buffer_ = NULL;
}

void GetURLHandler::Start() {
  pp::CompletionCallback cc =
      cc_factory_.NewCallback(&GetURLHandler::OnOpen);
  url_loader_.Open(url_request_, cc);
}

void GetURLHandler::OnOpen(int32_t result) {
  if (result != PP_OK) {
    ReportResultAndDie(url_, "pp::URLLoader::Open() failed", false);
    return;
  }
  
  // Check for Success 
  pp::URLResponseInfo response = url_loader_.GetResponseInfo();
  
  int32_t status_code = response.GetStatusCode();
  if( status_code < 200 || status_code > 299 ){
    std::stringstream ss;
    ss << status_code << std::endl;
    ReportResultAndDie(url_, "pp::URLLoader bad_http_status  ErrorCode:" +  ss.str(), false);
    return;
  }

  // Try to figure out how many bytes of data are going to be downloaded in
  // order to allocate memory for the response body in advance (this will
  // reduce heap traffic and also the amount of memory allocated).
  // It is not a problem if this fails, it just means that the
  // url_response_body_.insert() call in GetURLHandler::AppendDataBytes()
  // will allocate the memory later on.
  int64_t bytes_received = 0;
  int64_t total_bytes_to_be_received = 0;
  if (url_loader_.GetDownloadProgress(&bytes_received,
                                      &total_bytes_to_be_received)) {
    if (total_bytes_to_be_received > 0) {
      url_response_body_.reserve(total_bytes_to_be_received);
    }
  }
  // We will not use the download progress anymore, so just disable it.
  url_request_.SetRecordDownloadProgress(false);

  // Start streaming.
  ReadBody();
}

void GetURLHandler::AppendDataBytes(const char* buffer, int32_t num_bytes) {
  if (num_bytes <= 0)
    return;
  // Make sure we don't get a buffer overrun.
  num_bytes = std::min(READ_BUFFER_SIZE, num_bytes);
  // Note that we do *not* try to minimally increase the amount of allocated
  // memory here by calling url_response_body_.reserve().  Doing so causes a
  // lot of string reallocations that kills performance for large files.
  url_response_body_.insert(url_response_body_.end(),
                            buffer,
                            buffer + num_bytes);
}

void GetURLHandler::OnRead(int32_t result) {
  if (result == PP_OK) {
    // Streaming the file is complete, delete the read buffer since it is
    // no longer needed.
    delete [] buffer_;
    buffer_ = NULL;
    ReportResultAndDie(url_, url_response_body_, true);
  } else if (result > 0) {
    // The URLLoader just filled "result" number of bytes into our buffer.
    // Save them and perform another read.
    AppendDataBytes(buffer_, result);
    ReadBody();
  } else {
    // A read error occurred.
    ReportResultAndDie(url_,
                       "pp::URLLoader::ReadResponseBody() result<0",
                       false);
  }
}

void GetURLHandler::ReadBody() {
  // Note that you specifically want an "optional" callback here. This will
  // allow ReadBody() to return synchronously, ignoring your completion
  // callback, if data is available. For fast connections and large files,
  // reading as fast as we can will make a large performance difference
  // However, in the case of a synchronous return, we need to be sure to run
  // the callback we created since the loader won't do anything with it.
  pp::CompletionCallback cc =
      cc_factory_.NewOptionalCallback(&GetURLHandler::OnRead);
  int32_t result = PP_OK;
  do {
    result = url_loader_.ReadResponseBody(buffer_, READ_BUFFER_SIZE, cc);
    // Handle streaming data directly. Note that we *don't* want to call
    // OnRead here, since in the case of result > 0 it will schedule
    // another call to this function. If the network is very fast, we could
    // end up with a deeply recursive stack.
    if (result > 0) {
      AppendDataBytes(buffer_, result);
    }
  } while (result > 0);

  if (result != PP_OK_COMPLETIONPENDING) {
    // Either we reached the end of the stream (result == PP_OK) or there was
    // an error. We want OnRead to get called no matter what to handle
    // that case, whether the error is synchronous or asynchronous. If the
    // result code *is* COMPLETIONPENDING, our callback will be called
    // asynchronously.
    cc.Run(result);
  }
}

void GetURLHandler::ReportResultAndDie(const std::string& fname,
                                       const std::string& text,
                                       bool success) {
  ReportResult(fname, text, success);
  delete this;
}

void GetURLHandler::ReportResult(const std::string& fname,
                                 const std::string& text,
                                 bool success) {
  if (success){
    if (instance_) {
      pp::Var var_result("Successfully loaded: "+fname); 
      //instance_->PostMessage(var_result);
      file_provider_hook_->return_file_callback( text, false );    
    }
  }else{
    if (instance_) {
      pp::Var var_result("Error obtaining: "+fname);
      //instance_->PostMessage(var_result);
      file_provider_hook_->return_file_callback( text, true );    
    }
    
  }


}

